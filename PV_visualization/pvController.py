import subprocess

# write temporary args file for interfacing with the paraview-run sourceLoader

def setArgs(args):
    with open("args", 'w') as file:
        for arg in args:
            file.write(arg + ",")

# launching paraview gui

def launchPV(args):
    setArgs(args)
    subprocess.run(["paraview", "--script=./sourceLoader.py"])

# launching pvpython

def launchPVPython(args):
    setArgs(args)
    subprocess.run(["pvpython", "./sourceLoader.py"])

# testing

if __name__ == '__main__':
    launchPVPython(['test2.inp', 'test.npy'])
    #launchPV(['test2.inp', 'test.npy'])