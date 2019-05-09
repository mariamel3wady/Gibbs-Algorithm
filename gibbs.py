from random import *

def read_fasta(fp):
#FASTA file reader
#Do not change this function
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))
#============================================
# given file name retrive all Seq from file
filename = "data.fa"
allSeq = []
with open(filename) as fp:
    for name, seq in read_fasta(fp):
        allSeq.append(seq)
print("All data were read successfuly from: "+filename)
#============================================
#motifLen   given motif Len from the user
#seqLen     Seq. len from the FASTA file
#foundMotif List of all found motifs from all iterations
motifLen = int(input("Enter Desigred Motif Length."))
seqLen   = len(allSeq[0])
foundMotif = []
#============================================
#Main App
#
#
#Loop on all Given Seq
for i in range(len(allSeq)):
    currentSeq = allSeq[i] # Chosen Seq to cast away
    subSeq     = []        # Container of the reminder of the Sequnces
    subMotif   = []        # Container of all motifs taken randomly from Sequnces
    count      = []        # Count List to count the occurnce of every PROT/SEQ in all Sequnces
    prob       = []        # Probability of each Motif in Current iteration
    maxFrame = ""          # Frame with the highest probability
    maxFrameValue = -99999

    #------------------------------------------
    #Fill  subSeq list with all elements but the cast away
    for j in range(len(allSeq)):
        if j != i:
            subSeq.append(allSeq[j])
    #------------------------------------------
    #Fill subMotif List with previous motifs if found
    # Else with random motifs
    for j in range(len(foundMotif)):
        subMotif.append(foundMotif[j])
         
    for j in range(len(subSeq)- len(foundMotif)):
        rnd = randint(0,(seqLen-motifLen))
        subMotif.append(subSeq[j][rnd : rnd+motifLen])
    #------------------------------------------
    #Fill Count List with the count of every element in Col J
    for j in range(motifLen):
        tmp = {}
        for k in subMotif:
            value = k[j]
            try:
                tmp[value] = tmp[value] + 1
            except KeyError as e:
                tmp[value] = 1
        count.append(tmp)
    #Dividing the Counted value by Seq Len -1 to get the probability
    # of having this item in this Index
    # & adding one to remove all 0 values from the table
    for MAP in count:
        for key,val in MAP.items():
            MAP[key] = (val / (len(allSeq)-1))+1
    #------------------------------------------
    #Calulate the probability of each frame in cast away SEQ
    # & catching Max frame Probability
    for j in range(len(currentSeq)-motifLen):
        frame = currentSeq[j:j+motifLen]
        probV = 1
        for ch in frame:
            for k in count:
                try:
                   probV *= k[ch];
                except KeyError as e:
                    probV *= 1
        if probV > maxFrameValue:
            maxFrameValue = probV
            maxFrame      = frame
        prob.append(probV)
    foundMotif.append(maxFrame)
#============================================
#Print all Motifs found
print("Motifs : ")
print("\n".join(foundMotif))
#============================================
#Predicting the Best Motif by using Probabililty
count = []
for i in range(motifLen):
    tmp = {}
    for j in foundMotif:
        value = j[i]
        try:
             tmp[value] = tmp[value] + 1
        except KeyError as e:
            tmp[value] = 1
    for key,val in tmp.items():
        tmp[key] = (val / (len(allSeq)))
    count.append(tmp)

commonMotif = ""
for MAP in count:
    maxChar = ""
    maxCharValue = -999999
    for key,val in MAP.items():
        if val > maxCharValue:
            maxCharValue = val
            maxChar = key
    commonMotif += maxChar

print("Predicted Motif: ")
print(commonMotif)
            
        
    

              
        
    
    

    

        

        
        
            
        
    
    
            
    
    


