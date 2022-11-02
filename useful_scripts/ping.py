#! /usr/bin/env python3

import os
import sys
import time
import datetime
import matplotlib.pyplot as plt

#if this module is imported, one can simply run ping()

def loop_over_server(server):
  os.system(f'ping {server}.com -i 0.3 -c 100 >ping_{server}.txt')
  f=open(f"ping_{server}.txt",'r')
  g=open(f"statistics_{server}.txt",'w')
  string=''
  newstring=''
  numberx=''
  numbery=''
  numberz=''
  numberhangup=''
  X1=[]
  Y1=[]
  for lines in f:
      for letter in lines:
          if string=="time=" and letter!=" ":
              numbery+=letter
          elif string=="icmp_seq=" and letter!=" ":
              numberx+=letter
          elif letter==' ':
              string=''
          else:
              string+=letter
      if numberx !='' and numbery !='':
          X1.append(float(numberx))
          Y1.append(float(numbery))
          numberx=''
          numbery=''
      if "Request timeout for" in lines:
          for newletter in lines:
              if newletter !=' ':
                  numberhangup+=newletter
              else:
                  try:
                      X1.append(float(numberz))
                      Y1.append(None)
                  except ValueError:
                      pass
                  numberhangup=''
      if f"--- {server}.com ping statistics ---" in lines:
          g.writelines(lines)
          newstring='on'
      elif newstring=='on':
          g.writelines(lines)

  f.close()
  g.close()

  return X1, Y1

def ping():
  localtime=time.asctime( time.localtime(time.time() ))
  a= datetime.datetime.now()
  print(" ")
  print(" ---------------- ")
  print(" ")
  print("Program started to run at :", localtime)
                          
  X1,Y1 = loop_over_server("google")
  X2,Y2 = loop_over_server("gmail")
  X3,Y3 = loop_over_server("facebook")

  plt.plot(X1,Y1, marker='.', color='g', label='first 100 counts -> google.com')
  plt.plot(X2,Y2, marker='.', color='r', label='second 100 counts -> gmail.com')
  plt.plot(X3,Y3, marker='.', color='b', label='third 100 counts -> facebook.com')
  plt.title('Ping vs counts ')
  plt.xlabel('Counts')
  plt.ylabel('Ping (ms)')
  plt.legend(loc='best')
  plt.savefig('Ping_graph.pdf')
  plt.close()
  os.system("curl -L -o /dev/null ftp://speedtest.tele2.net/1000GB.zip 2>&1 |tr -u '\\r' '\\n' > bandwidth.log") # this link may be broken, another one has to be tested 
  time.sleep(60)
  os.system("killall 'curl'")

  os.system("cat bandwidth.log |awk '{print $NF}' >download_speed.log")
  #os.system("cat download_speed")
  #g=open("time.log")
  f=open("download_speed.log")
  
  data=[]
  t=[]
  timmm=0
  for i in f:
    #print "i=",i
    try:
      string=''
      mult=1./(1024*1024)
      for j in i:
        if j=="k":
          mult=1./1024
        elif j=="M":
          mult=1
        else:
          string=string+j
     # print "string=",string
      data.append(float(string)*mult)
      t.append(float(timmm))
      timmm+=1
    except ValueError:
      pass
  
  #for i in g:
  #    try:
  #        if "--" in i.split(":")[0]:
  #            t.append(0.0)
  #        else:
  #            t.append(float(i.split(":")[2])+60*float(i.split(":")[1])+3600*float(i.split(":")[0]))
  #    except:
  #        pass
      #print i
  
  #print "Time list:",t
  #print "Data list:",data
  
  plt.plot(t, data, marker='.', color='g', label='speed')
  plt.title('Speed vs time ')
  plt.xlabel('Time (s)')
  plt.ylabel('Download Speed (MB/s)')
  plt.legend(loc='best')
  plt.savefig('Bandwidth_graph.pdf')
  plt.close()
  datatemp=0
  for i in data:
    datatemp+=i
  timerun=t[-1]
  
  print("Average download speed (MB/s):", datatemp/timerun)
  #os.system("rm bandwidth.log bandwidth2.log download_speed.log time.log")
  
  
  b=datetime.datetime.now()
  localtime=time.asctime( time.localtime(time.time() ))
  c=b-a
  D=divmod(c.days*86400+c.seconds,60)
  print("Runtime :   ", D[0]," minutes and ", D[1], " seconds.")
  print("Program ended at :", localtime)
  print(" ")
  print(" ---------------- ")
  print(" ")

if __name__ == "__main__":
  ping()
