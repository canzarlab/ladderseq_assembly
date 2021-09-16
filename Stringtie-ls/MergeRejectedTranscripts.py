#!/usr/bin/env python
# coding: utf-8

# In[11]:


import HTSeq
import sys


# In[12]:


# rejectedGtf = "/home/schakraborty/Documents/HomeOffice/LadderSeq/refBasedAnalysis/30Mill/checkfilteredOut.gtf"
# acceptedGtf = "/home/schakraborty/Documents/HomeOffice/LadderSeq/refBasedAnalysis/30Mill/check.gtf"
# outputFileName = "/home/schakraborty/Documents/HomeOffice/LadderSeq/refBasedAnalysis/30Mill/transcripttobeputback.txt"

acceptedGtf = sys.argv[1]
rejectedGtf = sys.argv[2]
outputFileName = sys.argv[3]


# In[15]:


def getTranscriptIntronMapping(gtfFileName,threshold):
    gtfFile = HTSeq.GFF_Reader(gtfFileName, end_included=True)
    intronDict = dict()
    firstExonLength = -1
    lastExonLength = -1
    for feature in gtfFile:
        if feature.type == "transcript":
            if firstExonLength != -1 and lastExonLength != -1:
                if firstExonLength > threshold or lastExonLength > threshold :
                    intronDict[transcript_id] = intronList
                    
            transcript_id = feature.attr['transcript_id']
            start = 0
            end = 0
            intronList = list()
        elif feature.type == "exon":
            if(start == 0):
                start = feature.iv.end
                firstExonLength = abs(feature.iv.start - feature.iv.end)
                continue
            else:
                end = feature.iv.start
                intron = tuple((start,end))
                intronList.append(intron)
                start = feature.iv.end
                end = 0
                lastExonLength = abs(feature.iv.start - feature.iv.end)

    for key in list(intronDict.keys()):
        if(len(intronDict[key])==0):
            intronDict.pop(key)

    return(intronDict)


# In[16]:


def getIntronHashmap(gtfFileName):
    gtfFile = HTSeq.GFF_Reader(gtfFileName, end_included=True)
    intronHashmap = dict()
    for feature in gtfFile:
        if feature.type == "transcript":
            start = 0
            end = 0
        elif feature.type == "exon":
            if(start == 0):
                start = feature.iv.end
                continue
            else:
                end = feature.iv.start
                intron = tuple((start,end))
                if(intron not in intronHashmap):
                    intronHashmap[intron] = 1
                else:
                    intronHashmap[intron]+=1
                    
                start = feature.iv.end

    #print(intronHashmap)
    return(intronHashmap)


# In[17]:


rejectedIntronDict = getTranscriptIntronMapping(rejectedGtf,500)
acceptedIntronHashmap = getIntronHashmap(acceptedGtf)

transcriptBinaryDictionary = dict()
transcriptsToBePutBack = list()
#print(rejectedIntronDict)

for key,value in rejectedIntronDict.items():
    transList = list()
    count = 0
    novelIntron = False
    for v in value:
        if v not in acceptedIntronHashmap:
            transList.append(0)
            novelIntron = True
        else:
            count += 1
            transList.append(1)
    if novelIntron == True :
        transcriptsToBePutBack.append(key)
        
        
print(len(rejectedIntronDict))      
print(len(transcriptsToBePutBack))

with open(outputFileName, 'w') as filehandle:
    for listitem in transcriptsToBePutBack:
        filehandle.write('%s\n' % listitem)


# In[ ]:





# In[ ]:




