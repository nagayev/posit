from glob import glob
files=glob('lib/*.cpp')
files.extend(glob('lib/*.c'))
files.append('main.cpp')
print(files)
buf=""
for i in files:
    file=open(i)
    content=file.read()
    buf+="// "+file.name+"\n"
    buf+=content
with open('test.cpp','w') as output:
    output.write(buf)
