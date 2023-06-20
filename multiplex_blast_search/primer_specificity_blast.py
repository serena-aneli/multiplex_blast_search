import pandas
import numpy
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import csv
import pickle
import os.path
from pathlib import Path
import itertools

###TODO: SE VOGLIAMO METTERE AL POSTO DI 1 E 2 REDLABEL E REDNOT CONTROLLARE QUANDO USO INT

#################
### Read file ###
#################

def read_primer_file(file_path):
	primer_list = pandas.read_csv(file_path, names=['marker', 'type', 'seq'], sep="\t")
	return primer_list

####################
### Blast search ###
####################

# TODO: NON CREA LA CARTELLA DIRECTORY
def blast_search(directory, marker, type, primer_seq):
	directory = os.getcwd()
	file_name = Path("{directory}/{marker}_{type}.xml".format(directory = directory, marker = marker, type = type))

	if file_name.is_file():
				return "{directory}/{marker}_{type}.xml already downloaded".format(directory = directory, marker = marker, type = type)

	### blast search
	result_handle = NCBIWWW.qblast(program="blastn", database="nt", sequence=primer_seq, format_type='XML', short_query=True, megablast=True, hitlist_size = 500, expect = 1000) 

	### save results
	with open(file_name, 'w') as save_file: 
		blast_results = result_handle.read() 
		save_file.write(blast_results.replace('CREATE_VIEW', '')) #.replace('CREATE_VIEW', '') has been added to deal with the NCBI bug. The current version of Biopython on PyPy is not updated. When it will be updated, this part would not be necessary anymore.


##########################
### Alignments parsing ###
##########################
# funzione che parserizza e restituisce la tabella finale

def alignment_parsing(file_name):
	file_name = Path(file_name)
	with open(file_name) as saved_file:
		blast_records = list(NCBIXML.parse(saved_file))
	blast_records = list(blast_records)

	marker, type = os.path.splitext(file_name.name)[0].rsplit('_', maxsplit=1)

	### collect results and create table
	acc = []
	for record in blast_records:
		for alignment in record.alignments:
			for hsp in alignment.hsps:
				assert hsp.strand[0] == 'Plus'
				acc.append({
					'species': alignment.title.split("|")[4],
					'sbjct_id': alignment.hit_id,
					'sbjct_length': alignment.length,
					'aln_length': hsp.align_length,
					'aln_gaps': hsp.gaps,
					'aln_identities': hsp.identities,
					'aln_query_start': hsp.query_start,
					'aln_query_end': hsp.query_end,
					'aln_sbjct_start': hsp.sbjct_start,
					'aln_sbjct_end': hsp.sbjct_end,
					'aln_strand': hsp.strand[1],
				})
	tab = pandas.DataFrame.from_records(acc)
	tab[['gi', 'seq_name']] = tab['sbjct_id'].str.split("|", expand = True).iloc[:,[1,3]]
	file_excel = os.path.splitext(file_name.name)[0] + '.xlsx'	
	tab.to_excel(file_excel, index=False)
	#reds['red' + '_' + str(count) + '_' + primer] = tab


######################
### AMPLICON CHECK ###
######################

# al posto di fare rossi e unlabeled faccio tutti contro tutti

def find_amplicons(primers_file, output_file):
	
	primer_df = read_primer_file(primers_file)

	#output_file='/content/drive/MyDrive/FORENSIC_GENETICS/progetti/artefatti/amplicons_found_hitlist500.tsv'
	with open(output_file, 'w') as found_alignments: 

		print(
			'primer1', 'primer2', 'primer1_sequence', 'primer2_sequence', 'sequence_id', 'start', 'end', 'length', 'species', 
			'primer1_length', 'primer1_%_identity', 'primer1_gaps', 'primer1_start', 'primer1_end', 
			'primer2_length', 'primer2_%_identity', 'primer2_gaps', 'primer2_start', 'primer2_end', 'notes',
			sep="\t", file = found_alignments
		)

	### first iteration over primer pairs
		primer_df = primer_df.set_index(primer_df['marker'] + '_' + primer_df['type'].astype('str'))
		for i, j in itertools.combinations_with_replacement(primer_df.index, 2):
			print(i,j)
			df_i = pandas.read_excel(i + '.xlsx') 
			#if (i, j) != ('FGA_1', 'FGA_2'): continue
			df_j = pandas.read_excel(j + '.xlsx')
			### primer sequences
			i_seq = primer_df.at[i, 'seq']
			j_seq = primer_df.at[j, 'seq']

			common_findings = pandas.Series(list(set(df_i['seq_name']).intersection(set(df_j['seq_name']))))
			assert len(common_findings) == len(set(common_findings)) # check duplicates

			df_i_ = df_i.loc[df_i['seq_name'].isin(common_findings)]
			df_j_ = df_j.loc[df_j['seq_name'].isin(common_findings)]

			#print(len(common_findings))
			### second iteration over common alignments
			for al in common_findings:
				#if al != 'AF361104.2': continue
				df_i_al = df_i_.loc[df_i_['seq_name'] == al]
				df_j_al = df_j_.loc[df_j_['seq_name'] == al]

				### third iteration: each row of df_i_al with each row of df_j_al (if just an alignment is present over a certain sequence, it will do just one iteration)
				for index_i, row_i in df_i_al.iterrows():
					for index_j, row_j in df_j_al.iterrows():
						#print(row_i, row_j)
						if (row_i.at['aln_strand'] == 'Plus') and (row_j.at['aln_strand'] == 'Minus'):
							row_i.at['aln_sbjct_start']
							start = row_i.at['aln_sbjct_start']
							end = row_j.at['aln_sbjct_start']
							length = end - start

						elif (row_i.at['aln_strand'] == 'Minus') and (row_j.at['aln_strand'] == 'Plus'):
							start = row_j.at['aln_sbjct_start']
							end = row_i.at['aln_sbjct_start']
							length = end - start

						else: # plus-plus or minus-minus
							length = -1
						#from IPython import embed
						#embed()
						#assert False
						### check that the 3' of each primer is completely 'hybridized' (the query end must match with the length of the primer)
						#assert row_i['aln_query_end'] == len(i_seq)
						#assert row_j['aln_query_end'] == len(j_seq)
						if (row_i['aln_query_end'] != len(i_seq)) or (row_j['aln_query_end'] != len(j_seq)):
							notes = "mismatches at 3' of one primer"
						else:
							notes = ''
						#print(length, notes)
						if length > 0:
							#print(i, j, 'found:', al, start, end, length)
							print(
							i, j, i_seq, j_seq, al, start, end, length, row_i['species'],
							len(i_seq), row_i['aln_identities']/len(i_seq)*100, row_i['aln_gaps'], row_i['aln_query_start'], row_i['aln_query_end'],
							len(j_seq), row_j['aln_identities']/len(j_seq)*100, row_j['aln_gaps'], row_j['aln_query_start'], row_j['aln_query_end'], notes,
							sep="\t", file = found_alignments
							)
