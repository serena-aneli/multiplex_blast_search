parser = argparse.ArgumentParser(description='Primer blast search and amplicon size check')
parser.add_argument(
	'primer_list_path', metavar='primers', type=str, 
	help='path to the primer list files'
	)

args = parser.parse_args()


tab = read_primer_file(args.primer_list_path)

for index, row in tab.iterrows():
	print('searching {marker}, {pair} in BLAST'.format(marker=row['marker'], pair=row['type']))
	blast_search('blast_results', row['marker'], row['type'], row['seq']) # NON MI STAMPA IL MESSAGGIO CHE HA TROVATO
	search = row['marker'] + '_' + str(row['type'])
	alignment_parsing(search + '.xml')

find_amplicons('primers.txt', 'output_prova.tsv')
