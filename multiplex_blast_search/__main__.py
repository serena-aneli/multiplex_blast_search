import argparse

from .primer_specificity_blast import read_primer_file, blast_search, alignment_parsing, find_amplicons

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
		blast_search('blast_results', row['marker'], row['type'], row['seq']) # NON MI STAMPA IL MESSAGGIO CHE HA TROVATO
		search = row['marker'] + '_' + str(row['type'])
		alignment_parsing(search + '.xml')

	find_amplicons(primer_df, args.output_path)
