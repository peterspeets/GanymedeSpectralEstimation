import os
project_include_path  = os.path.dirname(os.path.realpath(__file__)) + "\\include\\"
denyList = ["BScan.h", "IO.h","IO.tpp","Settings.h","Tag.h","UtilityMathFunctions.h","UtilityMathFunctions.tpp","XML_Interpreter.h"]

for filename in os.listdir(project_include_path):
    
    if(filename in denyList):
        continue
    if(filename[:4] == "moc_" or filename[-2:] != '.h'):
        continue
    path = project_include_path + filename
    with open(path) as f:
        fileContents = f.read()
    
    if("Q_OBJECT" in fileContents):
        command ='C:\\Qt\\6.9.0\\mingw_64\\bin\\moc.exe {}{} -o {}moc_{}.cpp'.format(project_include_path,filename,project_include_path,filename[:-2])
        os.system(command) 
        
    