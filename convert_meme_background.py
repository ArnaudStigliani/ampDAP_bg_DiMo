#!/usr/bin/env python2.7

import sys


def convert_meme_background(methyl_bg, non_methyl_bg):	
	in_dict = {}
	comment_line = 1
	
	with open(methyl_bg, 'r') as in_bg:
		for line in in_bg:
			features = line.rstrip().split(' ', 1)
			if features[0] == "#":
				in_dict['comment{0}'.format(comment_line)] = features[1]
				comment_line += 1
			else:
				in_dict[features[0]] = float(features[1])
	
		
	with open(non_methyl_bg, 'w') as out_bg:
		out_bg.write("# {0}\n".format(in_dict['comment1']))
		out_bg.write("# {0}\n".format(in_dict['comment2']))
		out_bg.write("# {0}\n".format(in_dict['comment3']))
		out_bg.write("A {0:.3e}\n".format(in_dict['A']))
		out_bg.write("C {0:.3e}\n".format(in_dict['C']+in_dict['m']))
		out_bg.write("G {0:.3e}\n".format(in_dict['G']+in_dict['1']))
		out_bg.write("T {0:.3e}\n".format(in_dict['T']))
		out_bg.write("# {0}\n".format(in_dict['comment4']))
		out_bg.write("AA {0:.3e}\n".format(in_dict['AA']))
		out_bg.write("AC {0:.3e}\n".format(in_dict['AC']+in_dict['Am']))
		out_bg.write("AG {0:.3e}\n".format(in_dict['AG']+in_dict['A1']))
		out_bg.write("AT {0:.3e}\n".format(in_dict['AT']))
		out_bg.write("CA {0:.3e}\n".format(in_dict['CA']+in_dict['mA']))
		out_bg.write("CC {0:.3e}\n".format(in_dict['CC']+in_dict['mC']+in_dict['Cm']+in_dict['mm']))
		out_bg.write("CG {0:.3e}\n".format(in_dict['CG']+in_dict['mG']+in_dict['C1']+in_dict['m1']))
		out_bg.write("CT {0:.3e}\n".format(in_dict['CT']+in_dict['mT']))
		out_bg.write("GA {0:.3e}\n".format(in_dict['GA']+in_dict['1A']))
		out_bg.write("GC {0:.3e}\n".format(in_dict['GC']+in_dict['1C']+in_dict['Gm']+in_dict['1m']))
		out_bg.write("GG {0:.3e}\n".format(in_dict['GG']+in_dict['1G']+in_dict['G1']+in_dict['11']))
		out_bg.write("GT {0:.3e}\n".format(in_dict['GT']+in_dict['1T']))
		out_bg.write("TA {0:.3e}\n".format(in_dict['TA']))
		out_bg.write("TC {0:.3e}\n".format(in_dict['TC']+in_dict['Tm']))
		out_bg.write("TG {0:.3e}\n".format(in_dict['TG']+in_dict['T1']))
		out_bg.write("TT {0:.3e}\n".format(in_dict['TT']))
	


if __name__ == "__main__":
	
	methyl_bg = sys.argv[1]
	non_methyl_bg = sys.argv[2]
	convert_meme_background(methyl_bg, non_methyl_bg)
	

