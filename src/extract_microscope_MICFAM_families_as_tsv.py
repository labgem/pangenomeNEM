#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

import mysql.connector as mysql
import logging
import argparse
import wget
import os
import gffutils

logging.basicConfig(level = logging.DEBUG, format = '\n%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s', datefmt='%H:%M:%S')

if __name__=='__main__':

	parser = argparse.ArgumentParser(description='Extract Microscope families')
	parser.add_argument('-o', '--organisms', nargs='+', help="the O_id separated by spaces")
	parser.add_argument('-p', '--micfamparam', nargs=1, help="the MICFAM parameter (1=80id/80cov, 2=50id/80cov)")
	parser.add_argument('-f', '--families_file', nargs=1, type = argparse.FileType('w'), help="the out tsv file")
	parser.add_argument('-l', '--organism_list', nargs=1, type = argparse.FileType('w'), help="the out organism list")
	parser.add_argument('-g', '--gff_dir', nargs=1, type = str, help="the out direcotry containing gff files")

	options = parser.parse_args()

	if not os.path.exists(options.gff_dir[0]):
		os.makedirs(options.gff_dir[0])

	try:
	    conn = mysql.connect(host="mysqlagcdb.genoscope.cns.fr",port=3306,user="agc",passwd="admagc21",db="pkgdb_dev")  
	except:
	    print "Connexion with pkgdb_dev failed."
	    sys.exit(1)

	cur = conn.cursor()
	try:
		query = """CREATE TEMPORARY TABLE all_GO_id
	               SELECT GO_ori_id as GO_id
	               FROM Genomic_Object G
	               INNER JOIN Sequence S USING (S_id) 
	               INNER JOIN Replicon R USING (R_id)
	               INNER JOIN Organism O USING (O_id)
	               WHERE O_id IN (%s)
	               AND S_status = 'inProduction'
	               AND GO_type = 'CDS'
	               AND GO_update = 'current'
	               AND GO_status !=  'Artefact';""" % ",".join(options.organisms)
		logging.getLogger().debug(query)
		cur.execute(query)
	except mysql.Error as ex:
		logging.getLogger().error(ex)

	rows=[]

	try:
		query = """SELECT cluster_id, GO_id
	               FROM all_GO_id
	               INNER JOIN MICFAM_cluster USING(GO_id)
	               WHERE MICFAM_param_id=%s
	               ORDER BY cluster_id;""" % options.micfamparam[0]
		logging.getLogger().debug(query)
		cur.execute(query)
		rows = cur.fetchall()
	except mysql.Error as ex:
		logging.getLogger().error(ex)

	for row in rows:
		options.families_file[0].write(str(row[0])+"\t"+str(row[1])+"\n")

	PREFIX_URL="http://www.genoscope.cns.fr/agc/microscope/search/export.php?format=gff3&option=none&O_id="
	for org in options.organisms:
		url = PREFIX_URL+str(org)
		filename = wget.download(url, options.gff_dir[0]+"/"+str(org)+".gff")
		options.organism_list[0].write(str(org)+"\t"+str(filename))
		db_gff = gffutils.create_db(filename, ':memory:')
		for contig in db_gff.all_features(featuretype='region'):
			print(contig.seqid)
			print(contig.attributes['Note'])
			if contig.attributes['Note'].pop() == "chromosome circular complete sequence":
				options.organism_list[0].write("\t"+contig.seqid)
		options.organism_list[0].write("\n")
	conn.close()