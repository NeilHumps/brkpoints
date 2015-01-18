"""Name: Fragile_Finder.py
By: Neil Humphryes 18-Jan-2015
Purpose: To identify Fragile Site origins of Fused sequences resulting
from translocations
Input: Tab-delimited TICdb output
Output: de-duplicated csv file of most likely fragile sites"""

####################################################
# Modules
import optparse
import subprocess
import csv
from Bio.Blast import NCBIXML

####################################################
# Functions
def Parse_TICdbFile(TICdbFile):
    """Parse tab-delimited TICdb file and convert into fasta format
    because fasta is convenient blast input"""
    with open(TICdbFile) as file:
        contents = file.readlines()
    TIC_entries = [line.strip().split('\t') for line in contents]
    fasta = [ ]
    n=1
    for entry in TIC_entries :
        fasta.append(">" + '_'.join(entry[0:3]) + '\n' + entry[3])
    f=open('TIC.fasta' , 'w')
    f.write('\n'.join(i for i in fasta))
    f.close()

def BLAST_seqs(input_fastaFile, outfile_xml) :
    """Run Blast on fasta entries and output as XML, used blast to get best
        alignments rather than relying on perfect alignments/ allowing for
        mismatches using BWA"""
    # if running on a local server:
    subprocess.call(['blastn','-query',input_fastaFile,'-db',
                     '/usr/local/share/Blast/db/human_genomic', '-outfmt','5',
                     '-out', outfile_xml, '-num_alignments', '100'])

def Parse_BLAST(blast_XML_file) :
    """Parse BLAST XML output and compile relevant information"""
    BLAST_XML = open(blast_XML_file)
    blast_results = { }
    n=1
    # extract relevant data from XML
    for result in NCBIXML.parse(BLAST_XML) :
        name = str(n) + '_' + str(result.query)
        n=n+1
        length = result.query_length
        ALNs = [ ]
        for aln in result.alignments :
            if aln.hsps[0].expect < 5 :
                if aln.hit_def.find('Homo sapiens chromosome ') >= 0 :
                    ALNs.append([name, length, aln.hit_id, aln.hit_def,
                     aln.hsps[0].query_start, aln.hsps[0].query_end,
                     aln.hsps[0].sbjct_start, aln.hsps[0].sbjct_end,
                     aln.hsps[0].expect])
                else :
                    continue
        blast_results[name]=ALNs
    # compile best alignments into lists
    fragiles = [ ]
    for Alns in blast_results :
        best = { }
        # need to obtain best alignment for each side of fusion junction
        for Fragment in blast_results[Alns] :
            # Left-side alignment
            if Fragment[4] <= 2 :
                # only get best alignment, don't duplicate
                if best.get('left') :
                    continue
                else :
                    best['left'] = Fragment
                    # get Chromosome number
                    chr1 = str(Fragment[3]).split('Homo sapiens chromosome ')
                    chr2 = chr1[1].split(' ')
                    chr3 = chr2[0].split(',')
                    # determine strand and compile data
                    if Fragment[7] - Fragment[6] > 0 :
                        fragiles.append([Fragment[0]+'_a', chr3[0],
                                           Fragment[6], Fragment[7], 1])
                    else :
                        fragiles.append([Fragment[0]+'_a', chr3[0],
                                           Fragment[7], Fragment[6], -1])
            # Right-side alignment
            elif Fragment[5] >= Fragment[1]-5 :
                # only get best alignment, don't duplicate
                if best.get('right') :
                    continue
                else :
                    best['right'] = Fragment
                    # get Chromosome number
                    chr1 = str(Fragment[3]).split('Homo sapiens chromosome ')
                    chr2 = chr1[1].split(' ')
                    chr3 = chr2[0].split(',')
                    # determine strand and compile data
                    if Fragment[7] - Fragment[6] > 0 :
                        fragiles.append([Fragment[0]+'_b', chr3[0],
                                           Fragment[6], Fragment[7], 1])
                    else :
                        fragiles.append([Fragment[0]+'_b', chr3[0],
                                           Fragment[7], Fragment[6], -1])
    return fragiles

def remove_duplicates(fragile_list, distance) :
    """Remove duplicated sites, which are defined by sites that are located
        less than the specified distance in kb apart, default = 100kb"""
    x = int(distance)*1000
    by_chr = { }
    # organise fragile sites by chromosome
    for site in fragile_list[1:] :
        name = site[0].split('_')
        if name[4] == 'a' :
            Gene = name[1]
        elif name[4] == 'b' :
            Gene = name[2]
        if by_chr.get(Gene + '_' + site[1]) :
            by_chr[Gene + '_' + site[1]].append(site)
        else :
            by_chr[Gene + '_' + site[1]] = [site]
    # find duplicates
    dups = [ ]
    for chr in by_chr :
        positions = [ ]
        for i in by_chr[chr] :
            # make a list of sites within range of each site
            positions = [j for j in by_chr[chr] if
                     (int(i[2])+int(i[3]))/2 >= (int(j[2])+int(j[3]))/2 - x and
                      (int(i[2])+int(i[3]))/2 <= (int(j[2])+int(j[3]))/2 + x]
            # form one entry from each list covering full range of site
            dups.append([chr,
             '|'.join([positions[k][0] for k in range(0,len(positions))]),
             positions[0][1],
             min([min([int(positions[k][2]) for k in range(0,len(positions))]),
                  min([int(positions[k][3]) for k in range(0,len(positions))])
                  ]), max([
              max([int(positions[k][2]) for k in range(0,len(positions))]),
              max([int(positions[k][3]) for k in range(0,len(positions))])]),
             max([int(positions[k][4]) for k in range(0,len(positions))])])
    # remove duplicates
    no_dup = { }
    for site in dups :
        if no_dup.get(site[1]) :
            continue
        else :
            no_dup[site[1]] = site
    return no_dup

def write_csv(fragile_dict, out_File) :
    """Output fragile sites dictionary to csv file with header"""
    with open(out_File, 'w') as file:
        writer = csv.writer( file )
        writer.writerow(['Fragile_Site_Identifier', 'Fragile_site_entries',
                         'Chromosome','Start','End', 'Strand'])
        writer.writerows( [[fragile_dict[i][0],fragile_dict[i][1],
                            fragile_dict[i][2], fragile_dict[i][3],
                            fragile_dict[i][4], fragile_dict[i][5]]
                            for i in fragile_dict])
    
####################################################
# Main
if __name__ == '__main__' :
    # parse object for managing input options
    parser = optparse.OptionParser()

    # essential data, defines commandline options
    parser.add_option ('-t', dest = 'TICdbFile', default = '',
                       help = "This input is the tab-delimited \
                                TICdb output file.")
    parser.add_option ('-d', dest = 'distance', default = '100',
                       help = "This input is the threshold distance in kb \
                       to determine if two fragile sites are from the same \
                       fragile region.")
    parser.add_option ('-o', dest = 'out_File', default = 'Fragile_Sites.csv',
                       help = "This input is the name of the .csv \
                               output file.")
                   
    # load inputs
    (options, args) = parser.parse_args()

    # reads the inputs from commandline
    TICdbFile = options.TICdbFile
    out_File = options.out_File
    distance = options.distance
    
    # Parse input data
    Parse_TICdbFile(TICdbFile)
    # Run BLAST
    BLAST_seqs('TIC.fasta', 'TICdb_BLASTout.xml')
    # Parse BLAST output
    fragiles = Parse_BLAST('TICdb_BLASTout.xml')
    # remove duplicates
    no_dup = remove_duplicates(fragiles, distance)
    # write to .csv
    write_csv(no_dup, out_File)
