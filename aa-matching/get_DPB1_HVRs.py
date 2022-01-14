import json
import requests

# get hypervariable region array
url = "http://rest.hlatools.org/hladpb1/resources/hypervariableRegions"
response = requests.get(url)
# print(response.text)

hypervariableRegions = json.loads(response.text)

# pretty print JSON response with indenting
print(json.dumps(hypervariableRegions, indent = 4, sort_keys=True))


hvr_index = 0
for hvr in hypervariableRegions:
    # print (hvr)
    hvr_name = hypervariableRegions[hvr_index]['hypervariableRegionName']
    codon_list = hypervariableRegions[hvr_index]['codonNumberList']

    # get proteinSequenceList

    # get variantID

    hvr_index = hvr_index + 1
    print (hvr_name)
    print (codon_list)


# TODO - get list of alleles with each DPB1 HVR protein sequence motif
# some HVR variantIDs have more than one protein sequence
# assign variantID when it is defined
# there will be other AA motifs in IMGT/HLA without a variantID