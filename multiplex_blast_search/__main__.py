import sys
import argparse
from pathlib import Path

from .primer_specificity_blast import read_primer_file, blast_search, alignment_parsing, find_amplicons, NoAlignmentError

def main():
	parser = argparse.ArgumentParser(description='Primer blast search and amplicon size check')
	parser.add_argument(
		'primer_list_path', metavar='primers', type=str, 
		help='path to the primer list files'
		)
	parser.add_argument(
		'output_path', metavar='output', type=str, 
		help='name of the output file'
		)	

	args = parser.parse_args()

	primer_df = read_primer_file(args.primer_list_path)

	for index, row in primer_df.iterrows():
		print('searching {marker}, {pair} in BLAST'.format(marker=row['marker'], pair=row['type']))
		base_path = "{marker}_{type}".format(**row)
		#file_name = Path(f"{directory}/{marker}_{type}.xml".format(directory = directory, marker = marker, type = type))
		xml_path = Path(base_path + '.xml')

		try:
			if xml_path.is_file():
				print(xml_path, "already downloaded")
			else:
				blast_search(xml_path, row['seq'], True)
		except ConnectionError as e:
			print("Connection error while downloading results for primer {primer_name}. Let's try again.".format(primer_name=row['marker']), e)
			blast_search(xml_path, row['seq'], True)
		try:
			alignment_parsing(xml_path, base_path + '.xlsx')
		except NoAlignmentError:
			print("No alignments found, let's try without Megablast", file=sys.stderr)
			blast_search(xml_path, row['seq'], False)
			#sys.exit(1)

	find_amplicons(primer_df, args.output_path)
