#! /usr/bin/python
import sys
import string
import array
import struct
import re

# argument is name of data file (e.g. patient.001)
name = sys.argv[1]
n = open(name,mode='rb')

# read header
a = n.read(21)
print 'Version:',string.join(struct.unpack('21p',a),'')
a = n.read(21)
sampleRate = int(string.join(struct.unpack('21p',a),''))
print 'Sample rate:',sampleRate
a = n.read(21)
nbChannels = int(string.join(struct.unpack('21p',a),''))
print 'Nb of channels:',nbChannels
a = n.read(21)
experienceName = string.join(struct.unpack('21p',a))
print 'Experience Name:',experienceName
a = n.read(21)
print 'First Name:',string.join(struct.unpack('21p',a),'')
a = n.read(21)
print 'Sex:',string.join(struct.unpack('21p',a),'')
a = n.read(21)
print 'Birth date:',string.join(struct.unpack('21p',a),'')
a = n.read(21)
EPstudyNb = string.join(struct.unpack('21p',a),'')
EPstudyNb = 'NH'
print 'EP Study Nb:',EPstudyNb
a = n.read(21)
print 'Type of study:',string.join(struct.unpack('21p',a),'')
a = n.read(21)
fileName = string.join(struct.unpack('21p',a),'')
#fileName = 'PULLBACK_2_1635'
print 'Name of file including extension:',fileName
a = n.read(18)
print 'Time:',string.join(struct.unpack('18p',a),'')
a = n.read(4)
print '32bit system time in msec since start of program:',re.sub('\D','',`struct.unpack('i',a)`)
a = n.read(1)
print 'Gain of ECG channels:',re.sub('\D','',`struct.unpack('b',a)`)
a = n.read(22)
print 'VFG_info:',string.join(struct.unpack('22p',a),'')
a = n.read(21)
print 'Spare4:',string.join(struct.unpack('21p',a),'')
a = n.read(21)
print 'Spare5:',string.join(struct.unpack('21p',a),'')
a = n.read(21)
print 'Spare6:',string.join(struct.unpack('21p',a),'')
a = n.read(18)
print 'Spare7:',string.join(struct.unpack('18p',a),'')

# read channel labels and open files to write data
names = (['ECG', 'ECG','ECG','Pressure','Pressure'])
channelFiles = []
for i in range(128):
  a = n.read(8)
  if i < nbChannels:
    print 'Channel',i+1,'label:',string.join(struct.unpack('8p',a),'')
    #channelFiles.append(open(EPstudyNb+'-'+fileName+'_'+`i`+'.txt',mode='w'))
    channelFiles.append(open(names[i]+'_'+`i`+'.txt',mode='w'))

# go to data part
a = n.seek(8192)

#read data
eof = 0
samplenumber = 0

while not eof:
  samplenumber = samplenumber+1
  print 'Sample Number->',samplenumber
  for i in range(nbChannels):
    a = n.read(2)
    if a != "":
      print struct.unpack('2p',a)
      v = re.sub('\D','',`struct.unpack('H',a)`)
      channelFiles[i].write(v+'\n')
      print 'val:',v
      #print 'ChannelNo:',i,'NbChannels:',nbChannels
      
    else:
      eof = 1
#  if samplenumber == 45896:
#    eof=1
            
n.close()
for i in range(nbChannels):
  channelFiles[i].close()

