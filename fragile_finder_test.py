"""Name: fragile_finder_test.py
By: Neil Humphryes 18-Jan-2015
Purpose: To identify if fragile site alignments are actually proximal to the
documented gene loci
Input: de-duplicated csv file of most likely fragile sites
Output: % of correctly matched sites and
csv file with details of proximity to gene orf"""

####################################################
# Modules
import optparse
import subprocess
import csv
from Bio import SeqIO
from Bio import Entrez
Entrez.email = 'humphryesna@gmail.com'
from Bio.Blast import NCBIXML

####################################################
# Functions

def get_seqs(fragile_csv, range_start, range_end):
    """Obtain sequences for fragile site-associated genes from GenBank.
    range_start and range_end determine range of list entries to include,
    must be less than 100. Only 100bp of each gene sequence was obtained to
    make blast search faster - whole sequence not required"""
    # ensure entry range is compatible
    if int(range_end) - int(range_start) > 100 or \
            int(range_end) - int(range_start) < 1:
        range_start = 0
        range_end = 50
    # parse csv to get gene names
    with open(fragile_csv, 'r') as csvfile:
        fragiles_file = csv.reader(csvfile)
        fragiles = [row for row in fragiles_file]
    gene_seqs = []  
    for sites in fragiles[int(range_start):int(range_end)+1]:
        name = sites[0].split('_')
        gene = name[0]
        # retrieve gene info from Entrez database - get genbank identifier
        handle = Entrez.esearch(db='nucleotide',
                    term=str('Homo sapiens[Orgn] AND ' + gene + '[Gene]'))
        record = Entrez.read(handle)
        # use first record
        identifier = record['IdList'][0]
        handle = Entrez.efetch(db='nucleotide', id=identifier, rettype='gb',
                               retmode='text')
        record = SeqIO.read(handle, 'genbank')
        gene_seqs.append('>'+sites[1]+':'+str(len(record.seq))+'\n'+
                         str(record.seq[0:100]))
    f = open('Gene_seqs.fasta', 'w')
    f.write('\n'.join(i for i in gene_seqs))
    f.close()

def blast_seqs(input_fastafile, outfile_xml):
    """Run Blast on fasta entries and output as XML, used the same blast
    database to ensure genomic coordinates are compatible with previously-
    obtained coordinates. Only require best alignments so reduced number
    of alignments to 10"""
    # if running on a local server:
    subprocess.call(['blastn', '-query', input_fastafile, '-db',
                     '/usr/local/share/Blast/db/human_genomic', '-outfmt', '5',
                     '-out', outfile_xml, '-num_alignments', '10'])

def parse_blast_geneseq(blast_xml_file):
    """ Parse BLAST XML output and compile relevant information """
    blast_xml = open(blast_xml_file)
    blast_results = {}
    # extract relevant data from XML
    for result in NCBIXML.parse(blast_xml):
        alns = []
        for aln in result.alignments:
            if aln.hsps[0].expect < 5:
                if aln.hit_def.find('Homo sapiens chromosome ') >= 0:
                    alns.append([str(result.query), aln.hit_id, aln.hit_def,
                     aln.hsps[0].query_start, aln.hsps[0].query_end,
                     aln.hsps[0].sbjct_start, aln.hsps[0].sbjct_end,
                     aln.hsps[0].expect])
                else:
                    continue
        if len(alns) >= 1:
            blast_results[str(result.query)] = alns
    # compile best alignments into lists
    gene_pos = []
    for hits in blast_results:
        if len(blast_results[hits]) > 0:
            # get chromosome number
            chr1 = str(blast_results[hits][0][2]).split(
                'Homo sapiens chromosome ')
            chr2 = chr1[1].split(' ')
            chr3 = chr2[0].split(',')
            gene_pos.append([hits, chr3[0], blast_results[hits][0][5]])
        else:
            continue
    return gene_pos

def compare_loci(fragile_csv, gene_pos, distance, outfile_csv):
    """Compare position of genes mentioned in original TICdb file with fragile site position. Output statement displaying the % sites that are proximal to stated genes, and a csv file output showing distances between genes and fragile sites, and details of genes and sites that failed to colocalise according to previous threshold distance >100kb apart, and provide description"""
    with open(fragile_csv, 'r') as csvfile:
        fragiles_file = csv.reader(csvfile)
        fragiles = [row for row in fragiles_file]
    # make two dictionaries to organise gene and site data with
    # identical keys
    site_dict = {}
    for site in fragiles[1:]:
        site_dict[site[1]] = [site[2], (int(site[3])+int(site[4]))/2]
    gene_dict = {}
    for gene in gene_pos:
        header = gene[0].split(':')
        gene_dict[header[0]] = [gene[1], gene[2], header[1]]
    # make two lists for gene loci that match site loci, and one for
    # non-matches
    match = []
    different = []
    # determine if gene locus matches site locus based on threshold
    # distance.
    X = int(distance)*1000
    for gene in gene_dict:
        if site_dict.get(gene):
            if gene_dict[gene][0] == site_dict[gene][0]:
                if abs(gene_dict[gene][1] - site_dict[gene][1]) <= \
                       X + int(gene_dict[gene][2]):
                    match.append([gene, abs(gene_dict[gene][1] -
                                 site_dict[gene][1])])
                else:
                    different.append([gene, abs(gene_dict[gene][1] -
                                      site_dict[gene][1]),
                                      'Same chromosome but >100kb apart'])
            else:
                different.append([gene, 0, 'Different chromosomes'])
        else:
            different.append([gene, 0, 'Gene not in fragile site list'])
    # print statement showing % of matched sites
    # x = number of matched, y = number of different
    x = float(len(match))
    y = float(len(different))
    print str(int(x/(x+y)*100)) + \
     '% of the sites match their corresponding gene locus'
    # write lists to csv starting with non-matching sites first with comments
    with open(outfile_csv, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(['Fragile_Site_Identifier', 'Distance from Gene',
                         'Comments'])
        writer.writerows(different)
        writer.writerows(match)

####################################################
# Main
if __name__ == '__main__':
    # parse object for managing input options
    parser = optparse.OptionParser()

    # essential data, defines commandline options
    parser.add_option('-i', dest='fragile_csv',
                       default='Fragile_Sites.csv',
                       help="This input is the fragile site \
                               csv file.")
    parser.add_option('-s', dest='range_start', default='1',
                       help="This range specifies which fragile sites \
                               to test but must not exceed 100")
    parser.add_option('-e', dest='range_end', default='50',
                       help = "This range specifies which fragile sites \
                               to test but must not exceed 100")
    parser.add_option('-d', dest='distance', default='100',
                       help="This input is the threshold distance in kb \
                       to determine if the gene and fragile site colocalise.")
    parser.add_option('-o', dest='outfile_csv', default='fragile_test.csv',
                       help="This input is the name of the \
                               .csv output file.")
    # load inputs
    (options, args) = parser.parse_args()

    # reads the inputs from commandline
    fragile_csv = options.fragile_csv
    range_start = options.range_start
    range_end = options.range_end
    distance = options.distance
    outfile_csv = options.outfile_csv
    
    # Get gene sequences from GenBank
    get_seqs(fragile_csv, range_start, range_end)
    # Blast gene sequences to get genomic coordinates
    blast_seqs('Gene_seqs.fasta', 'Gene_seqs_BLASTout.xml')
    # Parse Blast output and get alignment positions
    gene_pos = parse_blast_geneseq ('Gene_seqs_BLASTout.xml')
    # compare gene positions to fragile site positions
    compare_loci (fragile_csv, gene_pos, distance, outfile_csv)



