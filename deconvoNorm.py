#!/usr/bin/python
# encoding: utf-8

# deconvoNorm
# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 1.0
# 2017-02-08
#
# Tested on : Ubuntu 14.04 LTS, 16.04 LTS
#      with : Python 2.7.10, 2.7.11


"""
deconvoNorm module
"""

import os
import sys
import argparse
from subprocess import call
import shutil
from math import ceil
from itertools import izip


def transition():
    """
    Transition between steps
    """
    print "-" * 60


##################################
#                                #
#####  Variable validation   #####
#                                #
##################################

def _validFile(path):
    """
    Intern function that raises an error if filename
    doesn't exists or is not in BAM format

    Return path
    """

    if not os.path.isfile(path):
        err = "The input file doesn't exist."
        raise argparse.ArgumentTypeError(err)
    else:
        try:
            open(path, 'r')
        except:
            err = "The input file can't be opened in reading mode."
            raise argparse.ArgumentTypeError(err)

        if not path.endswith(".bam"):
            err = "The input file must be in BAM format."
            raise argparse.ArgumentTypeError(err)
    return path


def _validDirectory(path):
    """
    Intern function that raises an error if filename
    doesn't exists or is not in narrowPeak format

    Return path
    """

    if not os.path.isdir(path):
        err = "The IP folder is not a directory."
        raise argparse.ArgumentTypeError(err)
    else:
        if not any(fname.endswith('.bam') for fname in os.listdir(path)):
            err = "The IP folder must contain bam files."
            raise argparse.ArgumentTypeError(err)
        else:
            return path


def _validFragmentLength(value):
    """
    Intern function that raises an error if
    the fragment length is below 1

    Return fragment length
    """

    try:
        value = int(value)
    except:
        raise argparse.ArgumentTypeError("The fragment length must be a int.")

    if value <= 0:
        raise argparse.ArgumentTypeError("The fragment length must be over 0.")

    return value


def _validWindowSize(value):
    """
    Intern function that raises an error if the
    window size is below 1 and an even number

    Return window size
    """

    try:
        value = int(value)
    except:
        raise argparse.ArgumentTypeError("The window size must be a int.")

    if value <= 0:
        raise argparse.ArgumentTypeError("The window size must be over 0.")

    if value % 2 == 0:
        raise argparse.ArgumentTypeError(
            "The window size must be an odd number.")

    return value


def _validThreshold(value):
    """
    Intern function that raises an error if the
    threshold value is below 0

    Return threshold
    """

    try:
        value = float(value)
    except:
        raise argparse.ArgumentTypeError("The threshold must be a float.")

    if value < 0:
        raise argparse.ArgumentTypeError("The threshold must be at least 0.")

    return value


def parseArgv(argv=None):
    """
    This function allows to parse the command line input options

    Return parsed input parameters
    """

    #######################
    ### Parse arguments ###
    #######################
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="  This script allows to use deconvolution"
                  + " to normalize sequencing data.",
        epilog="""
example 1 (one file):
  python deconvoNorm.py -i input.bam -f ip.bam -c MmrDNA -l 100 -w 25 -t 10 -o output

example 2 (many files):
  python deconvoNorm.py -i input.bam -d ipFolder -c MmrDNA -l 100 -w 25 -t 10 -o output
  """)

    ### OPTIONAL ARGUMENTS ###
    #parser.add_argument("--listchr", type=_validFile, metavar=('FILE'),
                        #help="List chromosomes from BAM file")

    ### FILES ###
    group1 = parser.add_argument_group("required files")

    # Input DNA file
    group1.add_argument("-i", "--input", type=_validFile, required=True,
                        help="Input DNA file (BAM format)")

    # IP file
    options = group1.add_mutually_exclusive_group(required=True)
    options.add_argument("-f", "--ipfile", type=_validFile,
                        help="Immunoprecipitation file (BAM format)")
    options.add_argument("-d", "--ipdirectory", type=_validDirectory,
                        help="Directory with many immunoprecipitation files")

    ### OPTIONS ###
    group2 = parser.add_argument_group("options")

    # Chromosome name
    group2.add_argument("-c", "--chrname", type=str, required=True,
                        help="Chromosome of interest (to normalize on)")

    # Fragment length
    group2.add_argument("-l", "--fragmentlength", default=100,
                        type=_validFragmentLength,
                        help="Sequenced fragment length (Default = 100)")

    # Window size
    group2.add_argument("-w", "--windowsize", default=25,
                        type=_validWindowSize,
                        help="Smoothing window size (Default = 25)")

    # Threshold
    group2.add_argument("-t", "--threshold", default=10,
                        type=_validThreshold,
                        help="Coverage threshold (Default = 10)")

    # Output name
    group2.add_argument("-o", "--outputname", type=str, required=True,
                        metavar=('OUTPUT'), help="Output folder name " +
                        "(if folder already exists, it will be overwritten)")

    # Keep intermediate files
    group2.add_argument("-k", "--keepfiles", action="store_true",
                        help="Keep intermediates files (Default = False)")

    # No RPM ratio
    group2.add_argument("-r", "--norpm", action="store_true",
                    help="Don't adjust in Read Per Million (Default = False)")

    return parser.parse_args(argv)


##################################
#                                #
##########  Functions   ##########
#                                #
##################################

def listChr(filename):
    """
    This function allows to list the chromosomes
    present in a BAM file

    Return the list of chromosomes
    """

    ### File validation ###
    if not os.path.isfile(filename):
        print("error : The input file doesn't exist.")
        sys.exit()
    try:
        open(filename, 'r')
    except IOError:
        print("error : The input file can't be opened in reading mode.")
        sys.exit()

    if not filename.endswith(".bam"):
        print("error : The input file must be in BAM format.")
        sys.exit()

    ### Create index if doesn't exist ###
    createIndex(filename)

    ### List chromosomes ###
    print "List of chromosomes for file {} :".format(filename)
    #out = check_output([])

    # Store idxstats in a temporary file
    with open("idxstats_deconvoNorm.txt", 'w') as f:
        call("samtools idxstats {}".format(filename), stdout=f, shell=True)

    # Read first column of idxstats
    with open("idxstats_deconvoNorm.txt", 'r') as r:
        for line in r:
            print(line.split()[0])

    os.remove("idxstats_deconvoNorm.txt")


def removeFiles(boolean, listOfFiles):
    """
    This function removes temporary
    files if -k option is OFF
    """

    ### If k is not present, removes the temporary files ###
    if not boolean:
        for oldFile in listOfFiles:
            os.remove(oldFile)


def createIndex(filename):
    """
    This function creates an index file
    if it doesn't already exists
    """

    basename = os.path.basename(filename)
    name = os.path.splitext(filename)[0]
    index = name + ".bai"

    # Check if file.bai or file.bam.bai exist
    if not os.path.isfile(index) and not os.path.isfile(filename + ".bai"):
        print "Creating index for file {}...".format(basename)
        call("samtools index {}".format(filename), shell=True)


def extractChr(filename, chrom, outFolder):
    """
    This function creates a bam file
    with only the chromosome of interest

    Return chromosome file names
    """

    ### Validation of chromosome name ###
    with open("idxstats_deconvoNorm.txt", 'w') as f:
        call("samtools idxstats {}".format(filename), stdout=f, shell=True)

    formattedChrom = "{}\t".format(chrom)
    if formattedChrom not in open("idxstats_deconvoNorm.txt", 'r').read():
        print("error : Chromosome {} not found in file {}"
              .format(chrom, filename))
        #print("Please, choose a chromosome in the list.")
        #print("To display the list of chromosomes names do :" +
        #      "python deconvoNorm.py --listChr filename.bam")
        sys.exit()

    os.remove("idxstats_deconvoNorm.txt")

    ### Creation of output name ###
    basename = os.path.basename(filename)
    newname = os.path.splitext(basename)[0] + "-{}.bam".format(chrom)
    outname = "{}/{}".format(outFolder, newname)

    ### Extraction of the chromosome of interest ###
    print("Extracting chromosome {} from file {}...".format(chrom, basename))
    call("samtools view -b {} {} > {}".format(filename, chrom, outname),
         shell=True)
    print("Temporary chromosome file created !")

    return outname


def bamtobed(filename):
    """
    This function creates a bed file
    with only the chromosome of interest

    Return bed file names
    """

    ### Creation of output name ###
    basename = os.path.basename(filename)
    outname = os.path.splitext(filename)[0] + ".bed"

    ### Bam to bed conversion ###
    print("Converting {} to bed format...".format(basename))
    call("bedtools bamtobed -i {} > {}".format(filename, outname), shell=True)
    print("Temporary bed file created !")

    return outname


def extReadLength(filename, fragmSize, chrLength):
    """
    This function extend read
    length from bed files

    Return extended bed file names
    """

    ### Creation of output name ###
    basename = os.path.basename(filename)
    outname = os.path.splitext(filename)[0] + "-l{}.bed".format(fragmSize)

    ### Extending read length ###
    print("Extending read length for file {}...".format(basename))

    # Read the bed file line by line
    with open(filename, 'r') as f, open(outname, 'w') as o:
        for line in f:
            # Assignation parameters
            chrom, start, end, name, score, strand = line.split()

            # Don't exceed the chromosome end
            if strand == "+":
                #newEnd = (chrLength) if (start + fragmSize > chrLength)
                #else (start + fragmSize)

                if (int(start) + fragmSize > int(chrLength)):
                    newEnd = chrLength
                else:
                    newEnd = int(start) + fragmSize

                o.write("{}\t{}\t{}\t{}\t{}\t{}\n"
                        .format(chrom, start, newEnd, name, score, strand))

            # Don't exceed the chromosome start
            else:
                #newStart = (0) if (end - fragmSize < 0) else (end - fragmSize)

                if (int(end) - fragmSize < 0):
                    newStart = 0
                else:
                    newStart = int(end) - fragmSize

                o.write("{}\t{}\t{}\t{}\t{}\t{}\n"
                        .format(chrom, newStart, end, name, score, strand))

    print("Temporary extended bed file created !")

    return outname


def genomecoverage(filename, chromInfo):
    """
    This function extract the
    coverage from bed files

    Return coverage bed file names
    """

    ### Creation of output name ###
    basename = os.path.basename(filename)
    outname = os.path.splitext(filename)[0] + "-cov.bed"

    ### Extending read length ###
    print("Extracting coverage for file {}...".format(basename))

    call("bedtools genomecov -i {} -g {} -d > {}"
         .format(filename, chromInfo, outname), shell=True)

    print("Temporary coverage bed file created !")

    return outname


def smooting(filename, winSize, threshold, chrLength):
    """
    This function smooths the coverage
    profile with a sliding window
    and replace by zero when the coverage
    falls below the threshold

    Return smoothed coverage bed file names
    """
    ### Creation of output name ###
    basename = os.path.basename(filename)
    outname = os.path.splitext(filename)[0] + "-w{}-t{}.bed".format(winSize, threshold)

    ### Smoothing coverage profile ###
    print("Smoothing coverage for file {}...".format(basename))

    slidingWindow = ceil(float(winSize) / 2)

    # Read the bed file line by line
    with open(filename, 'r') as f, open(outname, 'w') as o:
        lineNo = 0
        average = []
        for line in f:
            lineNo += 1
            # Assignation parameters
            chrom, pos, count = line.split()

            # Add the count to the end of the array
            average.append(float(count))

            # The first floor(winSize/2) lines will have a mean of zero
            if lineNo < slidingWindow:
                o.write("{}\t{}\t{}\n".format(chrom, pos, 0))

            # Once the average array contain winSize element, print the mean
            if len(average) == int(winSize):
                newPos = int(pos) - int(slidingWindow) + 1

                #newCount = sum(average) / len(average)

                #if newCount < threshold:
                #    newCount = 0

                mean = sum(average) / len(average)

                if float(mean) < float(threshold):
                    newCount = 0
                else:
                    newCount = mean

                o.write("{}\t{}\t{}\n".format(chrom, newPos, newCount))

                # Remove the first entry of the array
                del average[0]

            # The last floor(winSize/2) lines will have a mean of zero
            if lineNo == int(chrLength):
                newPos = int(pos) - int(slidingWindow) + 2
                for i in range(newPos, int(chrLength) + 1):
                    o.write("{}\t{}\t{}\n".format(chrom, i, 0))

    print("Temporary smoothed bed file created !")

    return outname


def scaleToRPM(filename, rpm):
    """
    This function scales the coverage
    in per million reads

    Return scaled coverage bed file names
    """

    ### Creation of output name ###
    basename = os.path.basename(filename)
    outname = os.path.splitext(filename)[0] + "-rpm.bed"

    print("Scaling to RPM for file {}...".format(basename))

    with open(filename, 'r') as f, open(outname, 'w') as o:
        for line in f:
            chrom, pos, count = line.split()

            newCount = float(count) * float(rpm)

            o.write("{}\t{}\t{}\n".format(chrom, pos, newCount))

    print("Temporary rpm bed file created !")

    return outname


def divisionOfInput(ipfile, inputfile):
    """
    This function divide the ip
    signal by the input signal

    Return division bed file names
    """

    basename = os.path.basename(ipfile)
    outname = os.path.splitext(ipfile)[0] + "_norm.bed"

    print("Division by input for file {}...".format(basename))

    # Scan the ip file and the input file simultaneously
    with open(ipfile, 'r') as ip, open(inputfile, 'r') as inp, open(outname, 'w') as o:
        for lineIp, lineInput in izip(ip, inp):
            # Get parameters
            chromIP, posIP, countIP = lineIp.split()
            chromInput, posInput, countInput = lineInput.split()

            # Avoid division by zero error
            if float(countInput) == 0.0:
                division = 0
            else:
                division = float(countIP) / float(countInput)

            o.write("{}\t{}\t{}\n".format(chromIP, posIP, division))

    print("Temporary divided bed file created !")

    return outname


def bedtobedgraph(filename):
    """
    This function convert bed file
    into bedgraph format
    """

    basename = os.path.basename(filename)
    outname = os.path.splitext(filename)[0] + ".bedgraph"

    print("Final conversion to bedgraph format for file {}...".format(basename))

    with open(filename, 'r') as f, open(outname, 'w') as o:

        # Header of bedgraph file
        o.write("track type=bedgraph visibility=full color=0,0,204 altColor=204,0,0 autoScale=on maxHeightPixels=128:128:11 viewLimits=0.0:25.0 yLineMark=0.0 windowingFunction=maximum\n")

        for line in f:
            chrom, end, count = line.split()

            start = int(end) - 1

            o.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, count))

    print("Final bedgraph file created !")


##################################
#                                #
##########  M A I N  #############
#                                #
##################################

def main():
    """
    Main function

    """

    ## List chromosomes ###
    #if sys.argv[1] == "--listChr":
        #listChr(sys.argv[2])
        #sys.exit
    #else:
        ## Parse argv
        #parseArgv()

    ### Parse argv ###
    argv = parseArgv()

    ### List immunoprecipitation files in an array ###
    listIPfiles = []
    if argv.ipdirectory:
        # Put the bam files in a list (except from input if present)
        for bamFile in os.listdir(argv.ipdirectory):
            path = "{}/{}".format(argv.ipdirectory, bamFile)
            if bamFile.endswith(".bam") and path != argv.input:
                listIPfiles.append(path)
        if len(listIPfiles) == 0:
            print("error : Only input present in {}".format(argv.ipdirectory))
            sys.exit()
    elif argv.ipfile:
        listIPfiles.append(argv.ipfile)

    ### New list with all the files ###
    listFiles = list(listIPfiles)
    listFiles.append(argv.input)

    ### Create index if necessary ###
    for bamFile in listFiles:
        createIndex(bamFile)

    ### Create output folder (overwrite if exists) ###
    out = argv.outputname
    if os.path.exists(out):
        shutil.rmtree(out)
    os.makedirs(out)

    ### Extract chromosome of interest ###
    transition()
    chrFiles = []
    for bamFile in listFiles:
        chrFiles.append(extractChr(bamFile, argv.chrname, out))

    ### Bam to bed conversion ###
    transition()
    bedFiles = []
    for chrFile in chrFiles:
        bedFiles.append(bamtobed(chrFile))

    # Remove temporary chromosomes files
    removeFiles(argv.keepfiles, chrFiles)

    ### Extension of read length ###
    # Extract chromosome length
    # Store idxstats in a temporary file
    with open("idxstats_deconvoNorm.txt", 'w') as f:
        call("samtools idxstats {}".format(argv.input), stdout=f, shell=True)

    # Find length of the chromosome of interest
    with open("idxstats_deconvoNorm.txt", 'r') as r:
        for line in r:
            chrom, length, aligned, notaligned = line.split()
            if chrom == argv.chrname:
                chrLength = length
                break

    # Extend read length
    transition()
    extFiles = []
    for bedFile in bedFiles:
        extFiles.append(extReadLength(bedFile, argv.fragmentlength, chrLength))

    # Remove temporary bed files
    removeFiles(argv.keepfiles, bedFiles)

    ### Extraction of coverage ###
    # Create chromInfo file
    with open("idxstats_deconvoNorm.txt", 'r') as r, open("chromInfo_deconvoNorm.txt", 'w') as o:
        for line in r:
            chrom, length, aligned, notaligned = line.split()
            if chrom == argv.chrname:
                o.write("{}\t{}\t{}\t{}"
                        .format(chrom, length, aligned, notaligned))
                break

    os.remove("idxstats_deconvoNorm.txt")

    # Extract coverage
    transition()
    covFiles = []
    for extFile in extFiles:
        covFiles.append(genomecoverage(extFile, "chromInfo_deconvoNorm.txt"))

    os.remove("chromInfo_deconvoNorm.txt")

    # Remove temporary extended bed files
    removeFiles(argv.keepfiles, extFiles)

    ### Smoothing coverage files ###
    transition()
    smoothFiles = []
    for covFile in covFiles:
        smoothFiles.append(smooting(covFile, argv.windowsize, argv.threshold, chrLength))

    # Remove temporary coverage bed files
    removeFiles(argv.keepfiles, covFiles)

    ### Scaling to RPM ###
    # Get RPM ratio associated with each file

    if not argv.norpm:
        dictRpm = {}
        for originalFile in listFiles:
            # Associate originalFile with extFile
            basename = os.path.basename(originalFile)
            new = os.path.splitext(basename)[0] + "-{}-l{}-cov-w{}-t{}.bed".format(argv.chrname, argv.fragmentlength, argv.windowsize, argv.threshold)
            outname = "{}/{}".format(argv.outputname, new)

            # Store idxstats in a temporary file
            with open("idxstats_{}.txt".format(basename), 'w') as f:
                call("samtools idxstats {}".format(originalFile),
                     stdout=f, shell=True)

            # Compute RPM ratio
            mappedReads = 0
            with open("idxstats_{}.txt".format(basename), 'r') as r:
                for line in r:
                    chrom, length, aligned, notaligned = line.split()
                    mappedReads += int(aligned)

                rpm = float(1000000) / float(mappedReads)

            # Add file and rpm to dictionnary
            dictRpm[outname] = rpm

            os.remove("idxstats_{}.txt".format(basename))

        # Scale to RPM
        transition()
        rpmFiles = []
        for smoothFile in smoothFiles:
            rpmFiles.append(scaleToRPM(smoothFile, dictRpm[smoothFile]))

        removeFiles(argv.keepfiles, smoothFiles)

    ### Division of IP files by input DNA file ###
    basename = os.path.basename(argv.input)

    if argv.norpm:
        new = os.path.splitext(basename)[0] + "-{}-l{}-cov-w{}-t{}.bed".format(argv.chrname, argv.fragmentlength, argv.windowsize, argv.threshold)
    else:
        new = os.path.splitext(basename)[0] + "-{}-l{}-cov-w{}-t{}-rpm.bed".format(argv.chrname, argv.fragmentlength, argv.windowsize, argv.threshold)

    newInput = "{}/{}".format(argv.outputname, new)

    transition()
    divFiles = []
    for ipfile in listIPfiles:

        basename = os.path.basename(ipfile)

        if argv.norpm:
            new = os.path.splitext(basename)[0] + "-{}-l{}-cov-w{}-t{}.bed".format(argv.chrname, argv.fragmentlength, argv.windowsize, argv.threshold)
        else:
            new = os.path.splitext(basename)[0] + "-{}-l{}-cov-w{}-t{}-rpm.bed".format(argv.chrname, argv.fragmentlength, argv.windowsize, argv.threshold)

        newIp = "{}/{}".format(argv.outputname, new)

        divFiles.append(divisionOfInput(newIp, newInput))

    # Remove temporary rpm or smooth bed files
    if argv.norpm:
        removeFiles(argv.keepfiles, smoothFiles)
    else:
        removeFiles(argv.keepfiles, rpmFiles)

    ### Final bedgraph file ###
    transition()
    for divFile in divFiles:
        bedtobedgraph(divFile)

    # Remove temporary coverage bed files
    removeFiles(argv.keepfiles, divFiles)

##################################

if __name__ == '__main__':

    main()

##################################
