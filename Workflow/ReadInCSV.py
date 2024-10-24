### this script is I/O for a output.csv file from REINVENT scaffold-decorator

def ReadCSV(PathToCSV):
    #### variables
    DecorationList=[]
    MoleculeList=[]
    Abundance=[]
    DecorDict=[]
    
    #### load data
    File=open(PathToCSV,'r')
    FileLines=File.readlines()
    File.close()
    
    for i in range(len(FileLines)-1):
        Line=FileLines[i+1].split(',')
        Abundance.append(int(Line[-1]))
        MoleculeList.append(Line[1])
        DistTmpList=[]
        if len(Line) == 5:
            if len(Line[3]) < 3:
                pass
            else:
                Deco=Line[3].split('\'')[1]
                if Deco in DecorationList:
                    DistTmpList.append(Line[3].split('\'')[1])
                    pass
                else:
                    DecorationList.append(Line[3].split('\'')[1])
                    DistTmpList.append(Line[3].split('\'')[1])
        if len(Line) == 6:
            Deco=Line[3].split('\'')[1]
            if Deco in DecorationList:
                DistTmpList.append(Line[3].split('\'')[1])
                pass
            else:
                DecorationList.append(Line[3].split('\'')[1])
                DistTmpList.append(Line[3].split('\'')[1])
            Deco=Line[4].split('\'')[1]
            if Deco in DecorationList:
                DistTmpList.append(Line[3].split('\'')[1])
                pass
            else:
                DecorationList.append(Line[4].split('\'')[1])
                DistTmpList.append(Line[3].split('\'')[1])
    
        DecorDict.append(DistTmpList)  
    return(DecorationList,MoleculeList,Abundance,DecorDict)
