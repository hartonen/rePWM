#!/usr/bin/env python

import argparse
import csv

from os import system

def rePWM():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    
    #MANDATORY PARAMETERS
    parser.add_argument("matrix",help="The input position frequency matrix file (tab-separeted).",type=str)
    parser.add_argument("seq",help="Full path to a Fasta-file containing sequences used to recalculate the PFM.",type=str,nargs='+')
    parser.add_argument("outfile",help="Full path to the output file (the recalculated PFM).",type=str)

    #OPTIONAL PARAMETERS, these are all parameters for MOODS

    group0 = parser.add_argument_group("Matrix file type:")
    group0.add_argument("--mtype",help="pwm or pfm (default=pfm).",type=str,choices=['pwm','pfm'],default='pfm')
    
    group1 = parser.add_argument_group("Threshold selection (exactly one required):")
    group1.add_argument("--p",help="Compute threshold from p-value.",type=float,default=None)
    group1.add_argument("--t",help="Use specified absolute threshold.",type=float,default=None)
    group1.add_argument("--B",help="Return at least the specified amount of best matches.",type=int,default=None)
    group1.add_argument("--Bforce",help="Return EXACTLY the specified amount of best matches to EACH of the masked motifs.",type=int,default=None)

    group2 = parser.add_argument_group("Search and model behaviour (optional):")
    group2.add_argument("--no-rc",help="If yes, disable matching versus the reverse complement strand (default=no).",type=str,choices=['yes','no'],default='no')
    group2.add_argument("--batch",help="If yes, do not recompute thresholds from p-values for each sequence separately (recommended when dealing with lots of short sequences, default=no).",type=str,choices=['yes','no'],default='no')
    group2.add_argument("--bg",help="Background distribution for computing thresholds from p-value with --batch (default is 0.25 for all alleles).",type=float,default=[0.25,0.25,0.25,0.25],nargs=4)
    group2.add_argument("--ps",help="Specify pseudocount for log-odds conversion (default=0.1).",type=float,default=0.1)
    group2.add_argument("--lo-bg",help="Background distribution for log-odds conversion (default is 0.25 for all alleles).",type=float,default=[0.25,0.25,0.25,0.25],nargs=4)

    parser.add_argument("--tmpdir",help="Directory used to store temporary files (masked PFMs/PWMs, MOODS output files), default = ./",type=str,default='./')
    parser.add_argument("--keeptmp",help="If yes, keep tmp files (default=no).",type=str,choices=['yes','no'],default='no')
    parser.add_argument("--v",help="If yes (=default), print progress to screen.",type=str,choices=['yes','no'],default='yes')
    args = parser.parse_args()

    if args.v=='yes': print "Creating masked matrices...",
    #first creating the masked versions of the input matrices
    matrix = []
    mlen = 0
    with open(args.matrix,'r') as csvfile:
        r = csv.reader(csvfile,delimiter='\t')
        for row in r:
            matrix.append([float(i) for i in row])
            mlen = len(row)

    #then saving each masked matrix into a temporary file
    for c in range(0,mlen):
        if args.mtype=='pwm':
            aux = 0.25
            auxfile = args.tmpdir+"rePWM_m"+str(c)+".pwm"
        else:
            aux = sum([row[c] for row in matrix])/4
            auxfile = args.tmpdir+"rePWM_m"+str(c)+".pfm"
        with open(auxfile,'w') as csvfile:
            w = csv.writer(csvfile,delimiter='\t')
            for row in matrix:
                if c==0: w.writerow([aux]+row[1:])
                elif c==len(row)-1: w.writerow(row[:-1]+[aux])
                else: w.writerow(row[:c]+[aux]+row[c:])
    if args.v=='yes': print "done!"

    #running MOODS for all the masked matrices
    moods_call = "moods_dna.py -o "+args.tmpdir+"moods.out "
    if args.mtype=='pwm':
        moods_call += "-S "
        for c in range(0,mlen): moods_call += args.tmpdir+"rePWM_m"+str(c)+".pwm "
    else:
        moods_call += "-m "
        for c in range(0,mlen): moods_call += args.tmpdir+"rePWM_m"+str(c)+".pfm "

    moods_call += "-s "
    for seq in args.seq: moods_call += seq+" "

    if args.p!=None: moods_call += "-p "+str(args.p)+" "
    elif args.t!=None: moods_call += "-t "+str(args.t)+" "
    elif args.B!=None: moods_call += "-B "+str(args.B)+" "
    #else: moods_call += "-B "+str(args.Bforce)+" "

    if args.no_rc=='yes': moods_call += "-R "
    if args.batch=='yes':
        moods_call += "--batch --bg "
        for l in args.bg: moods_call += str(l)+" "

    moods_call += "--ps "+str(args.ps)+" --lo-bg "
    for l in args.lo_bg: moods_call += str(l)+" "

    if args.v=='yes': print "running MOODS: "+moods_call+" ...",
    system(moods_call)
    if args.v=='yes': print "done!"

    #next we need to parse the results of individual masked matrices from the MOODS output file

    if args.v=='yes': print "Re-calculating the matrix...",
    new_matrix = []
    for i in range(0,4): new_matrix.append([0 for j in range(0,mlen)])

    #if exactly a set number of hits to each of the masked matrices is required, we remove all worse hits from moods.out
    if args.Bforce!=None:
        for c in range(0,mlen):
            if args.mtype=='pwm': name = "rePWM_m"+str(c)+".pwm"
            else: name = "rePWM_m"+str(c)+".pfm"

            #print 'grep "'+name+'" '+args.tmpdir+'moods.out | sort -t , -k5,5 -gr | head -'+str(args.Bforce)+' > '+args.tmpdir+'moods_sorted_'+str(c)+".out"

            print 'grep "'+name+'" '+args.tmpdir+'moods.out | sort -t , -k5,5 -gr > '+args.tmpdir+'moods_sorted_'+str(c)+"tmp.out"
            system('grep "'+name+'" '+args.tmpdir+'moods.out | sort -t , -k5,5 -gr > '+args.tmpdir+'moods_sorted_'+str(c)+"tmp.out")
            print 'head -'+str(args.Bforce)+' '+args.tmpdir+'moods_sorted_'+str(c)+"tmp.out"+' > '+args.tmpdir+'moods_sorted_'+str(c)+".out"
            system('head -'+str(args.Bforce)+' '+args.tmpdir+'moods_sorted_'+str(c)+"tmp.out"+' > '+args.tmpdir+'moods_sorted_'+str(c)+".out")
            #system('grep "'+name+'" '+args.tmpdir+'moods.out | sort -t , -k5,5 -gr | head -'+str(args.Bforce)+' > '+args.tmpdir+'moods_sorted_'+str(c)+".out")

        moodsfile = args.tmpdir+'moods_sorted.out'
        system('cat '+args.tmpdir+'moods_sorted_*.out > '+moodsfile)
    else: moodsfile = args.tmpdir+"moods.out"
            
    with open(moodsfile,'r') as csvfile:
        r = csv.reader(csvfile,delimiter=',')
        for row in r:
            seq = row[5]
            #ind = 0
            c = row[1].split('_')[1]
            c = int(c[1:-4])
            l = seq.upper()[c]
            if l=='A': new_matrix[0][c] += 1
            elif l=='C': new_matrix[1][c] += 1
            elif l=='G': new_matrix[2][c] += 1
            elif l=='T': new_matrix[3][c] += 1
            #ind += 1

    #saving the PFM into a file
    with open(args.outfile,'w') as csvfile:
        w = csv.writer(csvfile,delimiter='\t')
        for row in new_matrix: w.writerow(row)

    if args.v=='yes': print "done!"

    if args.keeptmp=='no':
        system("rm "+args.tmpdir+"moods.out")
        system("rm "+args.tmpdir+"rePWM_m*.p*m")
        if args.Bforce!=None:
            system("rm "+args.tmpdir+"moods_sorted_*.out")
            system("rm "+args.tmpdir+'moods_sorted.out')
        
                
#end

rePWM()
