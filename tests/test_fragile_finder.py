"""
Test of the fragile finder
"""

import pytest
import mock

import fragile_finder

@pytest.fixture
def input_file(request):
    """A mock input file of translocation sites"""
    return request.config.getoption("-t")

@pytest.fixture
def output_file(request):
    """A file like object I can write to and assert against"""
    return request.config.getoption("-o")

@pytest.fixture
def mock_blast_search(request):
    """A fake function we will call instead of the real blast query that
    will return known data"""
    print 'Mock BLAST search executed'

@pytest.mark.usefixtures('mock_blast_search')
def test_main(input_file, output_file):
    fragile_finder.main(input_file, output_file, 100)
    output_file.seek()
    results = output_file.readlines()
    assert len(results) > 0  # make more specific

def test_parse_blast():
    """test for correct function of fragile_finder.parse_blast by parsing 
    a test blast xml output to compare to the expected result"""
    # run parse_blast using test xml file
    fragiles_list = fragile_finder.parse_blast('fragile_finder_test_site_blastout.xml')
    # convert output to a single string
    fragiles_joined = [i[0]+'\t'+i[1]+'\t'+str(i[2])+'\t'+str(i[3])+'\t'+str(i[4]) 
                       for i in fragiles_list]
    fragiles = '\n'.join(i for i in fragiles_joined)
    # open expected output file as a single string
    with open('expected_parsed_blast.txt') as file:
        expected = file.read()
    # if strings are identical, parse_blast worked fine
    if fragiles == expected:
        print 'parse_blast worked correctly'
    else:
        print 'parsed BLAST does NOT match expected'
    #assert fragiles == expected

def test_remove_duplicates():
    """test for correct function of fragile_finder.remove_duplicates by 
    removing duplicates from a test file and comparing to test output"""
    with open('expected_parsed_blast.txt') as file:
        contents = file.readlines()
    parsed_blast = [line.strip().split('\t') for line in contents]
    nd = fragile_finder.remove_duplicates(parsed_blast, 100)
    # convert output to a single string
    no_dup_joined = [nd[i][0] + '\t' + nd[i][1] + '\t' + nd[i][2] + '\t' +
                     str(nd[i][3]) + '\t' + str(nd[i][4]) + '\t' + 
                     str(nd[i][5]) for i in nd]
    no_dup = '\n'.join(i for i in no_dup_joined)
    # open expected output file as a single string
    with open('without_duplicates.txt') as file:
        expected = file.read()
    # if strings are identical, parse_blast worked fine
    if no_dup == expected:
        print 'remove_duplicates worked correctly'
    else:
        print 'Duplicate removal does not match expected'
    
