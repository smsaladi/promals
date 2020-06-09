#!/usr1/local/bin/python

import os, sys, re, glob

# this script run promals given a fasta file
# process promals output alignment
# calculate conservation

name_width = 25

normal_finish1 = 0

include_res_num = 1

# read the log file
def read_promals_log(log_file, csv_string, caa_string, cutoff_csv):

	output_lines = ""

        # determine the maximum name length
        fp = open(log_file)
        maxlen = 0
        tmplen = 0
	record_begin = 0
        for line in fp:
                if re.match("CLUSTAL format", line): continue
		if re.match("ss:", line): record_begin = 1
		#if line.strip(): record_begin = 1
		if not record_begin: continue
		
		if "program finished" in line: 
			normal_finish = 1
			normal_finish1 = 1
			break
                list = line.split()
                if len(list) > 1: 
                        tmplen = len(list[0])
                        if tmplen > maxlen: maxlen = tmplen
        tmplen = len("Conservation:")
        if tmplen > maxlen: maxlen = tmplen

	record_begin = 0
	normal_finish = 0
	group_num = 0
        group_seq_num = 0
	ss_color_array1 = []
	ss_color_array2 = []
	ss_counts = 0
	fp = open(log_file)
	start_res_num = [0]
	end_res_num = [0]
	nblocks = 0
	seqnum = 0	
        blankline = 1
	
	for line in fp:
	
                if re.match("CLUSTAL format", line): continue
		if re.match("ss:", line): record_begin = 1
		#if line.strip(): record_begin = 1
		if not record_begin: continue
		
		if "program finished" in line: 
			normal_finish = 1
			normal_finish1 = 1
			break

                # determine the block number
                if line.strip():
                        if blankline == 1:
                                nblocks += 1
                                blankline = 0
                                # color and output conservation line
                                ss_array = line.split()[1]
                                ss_array_length = len(ss_array)
				if len(csv_string)==0: 
					#continue
					csv_string = '0' * ss_array_length
				this_csv = csv_string[0:ss_array_length]
                                #print this_csv
				csv_string = csv_string[ss_array_length:]
				this_color_csv = "Conservation:"
				#for i in range(14): this_color_csv += "&nbsp;"
				if include_res_num: 
					for i in range(maxlen-6): this_color_csv += " "
				else:
					for i in range(maxlen-12): this_color_csv += " "
				this_color_csv += color_csv_string(this_csv, cutoff_csv)
				output_lines += this_color_csv
				output_lines += "\n"
                else: blankline = 1
                        
		
		if re.match("ss:", line): 
			group_num += 1
			if(group_num%2==1): name_bg_color = 1
			else: name_bg_color = 0
				
			group_seq_num = 0
		
			#if group_num == 1: nblocks += 1
	
			ss_array = line.split()[1]
			ss_array_length = len(ss_array)
			ss_array = '-'+ss_array+'-';
			ss_color_array1 = []
			ss_color_array2 = []
			for i, c  in enumerate(ss_array):
				if(i==0): continue
				if(i==ss_array_length+1): continue
				# color helix red, and underscore
				if c=='H':
					if ss_array[i-1]!='H':
						#ss_color_array1.append("<FONT style='color:RED;'><u><b>")
						ss_color_array1.append("<FONT style='color:RED;'>")
					else: ss_color_array1.append("")
					if ss_array[i+1]!='H':
						#ss_color_array2.append("</b></u></FONT>")
						ss_color_array2.append("</FONT>")
					else: ss_color_array2.append("")
					continue
				# color strand blue, and underscore
				if c=='E':
					if ss_array[i-1]!='E':
						#ss_color_array1.append("<FONT style='color:BLUE;'><u><b>")
						ss_color_array1.append("<FONT style='color:BLUE;'>")
					else: ss_color_array1.append("")
					if ss_array[i+1]!='E':
						ss_color_array2.append("</FONT>")
						#ss_color_array2.append("</b></u></FONT>")
					else: ss_color_array2.append("")
					continue
				ss_color_array2.append("")
				ss_color_array1.append("")
			#for i in range(len(ss_color_array1)): print i, ss_color_array1[i], ss_color_array2[i]

			# color the csv string
			ss_counts += 1
                        '''
			if ss_counts == 1:
				if len(csv_string)==0: 
					#continue
					csv_string = '0' * ss_array_length
				this_csv = csv_string[0:ss_array_length]
				csv_string = csv_string[ss_array_length:]
				this_color_csv = "Conservation:"
				#for i in range(14): this_color_csv += "&nbsp;"
				if include_res_num: 
					for i in range(maxlen-6): this_color_csv += " "
				else:
					for i in range(maxlen-12): this_color_csv += " "
				this_color_csv += color_csv_string(this_csv, cutoff_csv)
				output_lines += this_color_csv
				output_lines += "\n"
                                '''
	
		elif line!='\n':
			name, seq = line.split()

			group_seq_num += 1
			seqnum += 1
			if(nblocks == 1) and ("Consensus_ss:" not in name):
                                fragment_length = len(seq.replace("-", ""))
                                if fragment_length>0: start_res_num.append(1)
                                else: start_res_num.append(0)
                                end_res_num.append(fragment_length);
				#start_res_num.append(1)
				#end_res_num.append(len(seq.replace("-", ""))-1)
				#print seqnum, nblocks, start_res_num[seqnum], end_res_num[seqnum]
			elif "Consensus_ss:" not in name:
                                #print line,"|"
                                fragment_length = len(seq.replace("-", ""))
                                if fragment_length > 0:
                                        start_res_num[seqnum] = end_res_num[seqnum]+1
                                else:
                                        start_res_num[seqnum] = end_res_num[seqnum]
				#start_res_num[seqnum] = end_res_num[seqnum]+1
				end_res_num[seqnum] = end_res_num[seqnum]+len(seq.replace("-", ""))
				#print seqnum, nblocks, start_res_num[seqnum], end_res_num[seqnum]
	
			#name = name[0:name_width]
                        for i in range(maxlen - len(name)): name = name + " "
			#name = "%-25s" %name
			#name = re.sub(" ", "&nbsp;", name)
			#if name_bg_color==0: name = '<span style="background-color: RGB(175,175, 175)">'+name+'</span>'
			#else: name = '<span style="background-color: RGB(235,235,235)">'+name+'</span>'
			# this is for representative
			if(group_seq_num==1) and ss_color_array1 and ss_color_array2:
				# make the name in bold and color in magenta
                                #name = '<b><span style="color: #800000">'+name+'</span></b>'
                                #name = '<span style="color: #800000">'+name+'</span>'
                                name = '<span style="color: #FF00FF">'+name+'</span>'
				seq1 = seq
				seq = ""
                                #print name, seq1
				for i, c in enumerate(seq1):
					seq += ss_color_array1[i] + seq1[i] + ss_color_array2[i]
			#name += "&nbsp;&nbsp;&nbsp;"
			#name += "   "

			# for consensus secondary structure
			if("Consensus_ss:" in name):
				#seq = re.sub("e", '<span style="color:BLUE">e</span>', seq)
				#seq = re.sub("h", '<span style="color:RED">h</span>', seq)
				seq = re.sub("\.", " ", seq)
	
			#web_line = name+seq+"<br>"
			if(include_res_num) and ("Consensus_ss:" not in name):
                                fragment_length = len(seq.replace("-", ""))
                                if fragment_length == 0:
                                        web_line = "%s %4s  %s %4s" %(name, ' ', seq, ' ');
                                else: web_line = "%s %4d  %s %4d" %(name, start_res_num[seqnum], seq, end_res_num[seqnum]);
				#web_line = "%s %4d  %s %4d" %(name, start_res_num[seqnum], seq, end_res_num[seqnum]);
                        elif include_res_num and ("Consensus_ss:" in name):
                                if(caa_string):
                                        name_aa = name_aa1 = "Consensus_aa:"
                                        name_aa1 = "<a href=http://prodata.swmed.edu/promals3d/info/consensus.html>" + name_aa1 + "</a>"
                                        for i in range(maxlen - len(name_aa)): name_aa1 = name_aa1 + " "
                                        name_aa = name_aa1
                                        colored_caa = color_consensus(caa_string[0:len(seq)])
                                        web_line = (name_aa+"       "+colored_caa)
                                        caa_string = caa_string[len(seq):]
                                        web_line += "\n"
                                name = re.sub("Consensus_ss:", "<a href=http://prodata.swmed.edu/promals3d/info/consensus_ss.html>Consensus_ss:</a>", name)
                                web_line += (name+"       "+seq)
                        elif not include_res_num and ("Consensus_ss:" in name): 
                                if(caa_string):
                                        name_aa = name_aa1 = "Consensus_aa:"
                                        name_aa1 = "<a href=http://prodata.swmed.edu/promals3d/info/consensus.html>" + name_aa1 + "</a>"
                                        for i in range(maxlen - len(name_aa)): name_aa1 = name_aa1 + " "
                                        name_aa = name_aa1
                                        colored_caa = color_consensus(caa_string[0:len(seq)])
                                        web_line = (name_aa+" "+colored_caa)
                                        caa_string = caa_string[len(seq):]
                                        web_line += "\n"
                                name = re.sub("Consensus_ss:", "<a href=http://prodata.swmed.edu/promals3d/info/consensus_ss.html>Consensus_ss:</a>", name)
                                web_line += (name+" "+seq)
                        else:
                                web_line = (name+" "+seq)
	
			output_lines +=  web_line
			output_lines += "\n"
	
		else:
			#output_lines += "<br>"
			output_lines += "\n"
			group_num = 0
			ss_counts = 0
			seqnum = 0

	fp.close()

	if record_begin: return output_lines
	fp = open(log_file)
	seqnum = 0
	for line in fp:
	
		if re.match("nblocks:", line): 
			record_begin = 1
			continue
		if not record_begin: continue
		
		if re.match("------------------", line): 
			normal_finish = 1
			normal_finish1 = 1
			break
		if "program finished" in line: 
			normal_finish = 1
			normal_finish1 = 1
			break

		#output_lines += line.replace(" ", "&nbsp;")
		if line!= '\n':
			seqnum += 1
		else:
			seqnum = 0
		if seqnum == 1:
			ss_array_length = len(line.split()[1])
			if len(csv_string)!=0:
				this_csv = csv_string[0:ss_array_length]
				csv_string = csv_string[ss_array_length:]
				this_color_csv = "Conservation:"
				#for i in range(14): this_color_csv += "&nbsp;"
				if include_res_num:
					for i in range(maxlen-7): this_color_csv += " "
				else:
					for i in range(maxlen-12): this_color_csv += " "
				this_color_csv += color_csv_string(this_csv, cutoff_csv)
				output_lines += this_color_csv
				output_lines += "\n"
			else:
				output_lines += "Conservation: none\n"

		if seqnum >= 1:
			name, seq = line.split()
                        for i in range(maxlen - len(name)): name = name + " "
			#name = "%-25s" %(name[0:25])
			#line = '<span style="background-color: RGB(175,175, 175)">%s</span>   %s\n' %(name, seq)
			line = '%s   %s\n' %(name, seq)
			
		output_lines += line
		#output_lines += "<br>"

	fp.close()
	return output_lines
		

def color_consensus(caa_string):

        color_str = ""
        for i in caa_string:
                if i == '.': color_str += '.'
                #elif i.isupper(): color_str = color_str + "<b><span style=\"background: black; color: white\">" + i + "</span></b>"
                elif i.isupper():
                        if i in "": #'KRH':
                                color_str = color_str + "<b><span style=\"color: blue\">" + i + "</span></b>"
                        elif i in "": #'DE':
                                color_str = color_str + "<b><span style=\"color: red\">" + i + "</span></b>"
                        else: color_str = color_str + "<b>" + i + "</b>"
                #elif i == 'h': color_str = color_str + "<span style=\"background: yellow\">"+i+"</span>"
                elif i in 'hl': color_str = color_str + "<span style=\"color: #347C17\"><i>"+i+"</i></span>"
                elif i == '+': color_str = color_str + "<b><span style=\"color: blue\">"+i+"</span></b>"
                elif i == '-': color_str = color_str + "<b><span style=\"color: red\">"+i+"</span></b>"
                #elif i == 'a': color_str = color_str + "<span style=\"color: #F87431\">"+'@'+"</span>"
                elif i == 'a': color_str = color_str + "<span style=\"color: #347C17\"><i>"+'@'+"</i></span>"
                #elif i in 'ts': color_str = color_str + "<span style=\"color: #F87431\">"+i+"</span>"
                else: color_str = color_str +  i 
        #print color_str
        return color_str

def color_csv_string(conservation_array, cutoff):

	span_bg_color = []
	j = 0

	blue_spectrum = [ "EEFFFF", "C6EFF7", "94D6E7", "63C6DE", "31B5D6", "00A5C6", "0084A5", "006B84", "005263", "00394A"]

	for i in range(0, 10):

		if i<cutoff: color_string = ' '
		elif i <=7: 
                        color_string = '<span style="color: #%s">%d</span>' %("000000", i)
		elif i >= 8:
                        color_string = '<span style="color: #%s">%d</span>' %("800000", i)	
		span_bg_color.append(color_string)
		continue
			
		if i<=5: color_string = " "
                #if i>=5: color_string = '<span style="color: #%s">%d</span>' %(blue_spectrum[i], i)
                if i>5 and i<=7: color_string = '<span style="color: #%s">%d</span>' %("000000", i)
                if i>=8: color_string = '<span style="color: #%s">%d</span>' %("800000", i)
		#if i>=8: color_string = '<span style="background-color: #99FF33">' + color_string + '</span>'
		#if i>=8: color_string = '<span style="background-color: #66FF33">' + color_string + '</span>'
		span_bg_color.append(color_string)
		continue

		if(i<=5): color_string = '<span style="background-color: #%s"> </span>' %(blue_spectrum[i])
		else: color_string = '<span style="background-color: #%s">%d</span>' %(blue_spectrum[i], i)
		#if(i==6 or i==7): color_string = '<b>'+color_string+'</b>'
                #if(i==8): color_string = '<b><span style="color: #BB00BB">'+color_string+'</span></b>'
                #if(i==9): color_string = '<b><span style="color: #FF0000">'+color_string+'</span></b>'

		#if(i==6 or i==7): color_string = color_string
                if(i==8): color_string = '<span style="color: #FFFFFF">'+color_string+'</span>'
                if(i==9): color_string = '<span style="color: #FFFFFF">'+color_string+'</span>'
                #if(i>=5): color_string = '<span style="color: #FFFFFF">'+color_string+'</span>'
		span_bg_color.append(color_string)
	#for i in span_bg_color: print i

	web_conservation_array = ""
	for i in conservation_array:
		int_value = int(i)
		web_conservation_array += span_bg_color[int_value]
                #print int_value, span_bg_color[int_value]
	#web_conservation_array += "<br>"
        #print web_conservation_array
	return web_conservation_array

def get_csv_string(alnfile, caa_freq):
	#print alnfile
        csv_string = ""
        csvfile = alnfile + ".csv"
        #print "caa_freq: ", caa_freq
        if os.path.isfile(csvfile):
                fp = open(csvfile)
                csv_string = fp.readline()
                caa_string = fp.readline()
                caa_level_string = fp.readline()
                #print csvfile
                #print caa_level_string
                caa_level_string = re.sub("consensus_level:\s+", "", caa_level_string)
                #print "consensus level: ", caa_level_string
                if csv_string: 
                        try: 
                                a1 = float(caa_level_string)
                                if a1 == caa_freq: return csv_string, caa_string
                        except:
                                pass

	command = "/home/jpei/t2/bar/promals_package/bin/al2co_consensus -i %s -t %s.csv.aln -g 0.25 -b 1000000 -consensus T -freq_cutoff %f > /dev/null" %(alnfile, alnfile, caa_freq)
        #print command
	#command = "/home/jpei/t2/bar/promals_package/bin/al2co -i %s -t /tmp/csv.aln -g 0.25 -b 1000000 > /dev/null" %(alnfile)
	os.system(command)
        if not os.path.isfile(alnfile+".csv.aln"):
        #if not os.path.isfile("/tmp/csv.aln"):
		return ""
        conservation_line = ""
        consensus_aa_line = ""
	for line in open(alnfile+".csv.aln").readlines():
	#for line in open("/tmp/csv.aln").readlines():
		if(re.match("Conservation:", line) ) : conservation_line = line
                if(re.match("Consensus_aa:", line) ) : consensus_aa_line = line
	csv_string = conservation_line.split()[1]
        caa_string = consensus_aa_line.split()[1]
        os.system("rm -f " + alnfile+".csv.aln")
        #print "consensus not match, recalculate consensus"
	return csv_string, caa_string

def run_promals(inputfile):
	
	command = "/home/jpei/profile_hmm/profile_hmm/progress %s -outfile %s.promals.aln > %s.log" %(inputfile, inputfile, inputfile)
	os.system(command)

	if(not os.path.isfile(inputfile+".promals.aln") ):
		print "No alignment generated"
		sys.exit(0)

if __name__ == "__main__":

	inputfile = sys.argv[1]
	#run_promals(inputfile)	
	#os.system("rm -f " + inputfile + ".csv.aln")
	out_log_file = inputfile.replace(".result", "")
	out_log_file = re.sub("\/php", "/QUERY_", out_log_file) + ".lock"
	if(len(sys.argv) > 2): 
		if sys.argv[2][0]!='-': out_log_file = sys.argv[2]
	cutoff_value = 5
        alignedorder = 1
        caa_freq = 0.8
	for i in range(len(sys.argv)):
		if sys.argv[i] == '-cutoff':
			cutoff_value = int(sys.argv[i+1])
		if sys.argv[i] == '-caa_freq':
			caa_freq = float(sys.argv[i+1])
		if sys.argv[i] == '-resnum':
			include_res_num = int(sys.argv[i+1]) # 1 or 0
		if sys.argv[i] == '-inputorder':
			alignedorder = 0
		
	mycsv, mycaa = get_csv_string(inputfile, caa_freq)
	if cutoff_value >=10: cutoff_value = 9
	#out_log_file = re.sub("\/php", "/QUERY_", out_log_file) + ".err"
	check_lock_finish = open(out_log_file)
	check_finish = 0
	for myline in check_lock_finish:
		if "program finished" in myline:
			check_finish = 1
        check_lock_finish.close()
	#if check_finish==0: sys.exit()
	web_alignment = read_promals_log(out_log_file, mycsv, mycaa, cutoff_value)
	#out_web_file = "/tmp/" + inputfile+".html"
	out_web_file = inputfile+".html"
	out_web_fp = open(out_web_file, "w")

	out_web_fp.write('<html><head><title>PROMALS3D Result</title></head><body>')
        out_web_fp.write('<h3 style="font-family: verdana"> <a href="http://prodata.swmed.edu/promals3d/info/promals3d_help.html#output1">Colored</a> PROMALS3D alignment</h3><hr>\n<pre style="font-family: \'Courier New\'">')
        '''
        if alignedorder:
                out_web_fp.write('<h3 style="font-family: verdana"> <a href="http://prodata.swmed.edu/promals3d/info/promals_output.html">Colored</a> PROMALS3D alignment (sequences in aligned order)</h3><hr><pre>')
        else:
                out_web_fp.write('<h3 style="font-family: verdana"> <a href="http://prodata.swmed.edu/promals3d/info/promals_output.html">Colored</a> PROMALS3D alignment (sequences in input order)</h3><hr><pre>')
        '''

	#############################
	'''
	ref_file = "/home/jpei/prefab4/ref/" + inputfile
	reflines = open(ref_file).readlines()
	for line in reflines: 
		tmpline = ""
		for aletter in line:
			if aletter.isupper():
				tmpline = tmpline+"<b><u>"+aletter+"</u></b>"
			else: tmpline += aletter
		tmpline += "<br>"
		out_web_fp.write(tmpline)
	out_web_fp.write("<br>")
	'''
	#############################
	
	out_web_fp.write(web_alignment)
	out_web_fp.write('<br><br>\n')
	#out_web_fp.write('<a href="http://prodata.swmed.edu//promals/info/promals_output.html">Here is information about this alignment format</a>')
	out_web_fp.write("</body></html>\n")
	out_web_fp.write("<!---promals_now_finished--->\n")

	out_web_fp.close()
	#print web_alignment

	'''
	if not normal_finish1:
		command = "rm -f %s" %out_web_file
		os.system(command)
		print "Here:"
	'''

				
