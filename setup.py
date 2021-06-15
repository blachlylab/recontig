from pyd.support import setup, Extension
import distutils.ccompiler
from distutils.command.clean import clean
import subprocess
import glob
import sys
import os
import shutil

projName = 'recontig'

htslibVersion = "1.12"

# default htslib path is local
htslibLocalPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),"htslib-{}".format(htslibVersion))

SHARED_EXT = ".so"
if sys.platform == "darwin":
    SHARED_EXT = ".dylib"
# download htslib from github
def getHtslib():
    sys.stdout.flush()

    try:
        retcode = subprocess.call(
            ["curl", "https://github.com/samtools/htslib/releases/download/{v}/htslib-{v}.tar.bz2".format(v=htslibVersion),"-o","htslib-{v}.tar.bz2".format(v=htslibVersion)]
            )
        if retcode != 0:
            return False
        else:
            return True
    except OSError as e:
        return False

# untar and decompress htslib package
def unpackHtslib():
    sys.stdout.flush()

    try:
        retcode = subprocess.call(
            ["tar", "-xjf","htslib-{}.tar.bz2".format(htslibVersion)],
            )
        if retcode != 0:
            return False
        else:
            return True
    except OSError as e:
        return False

# run configure
def runConfigure(option):
    sys.stdout.flush()

    try:
        retcode = subprocess.call(
            " ".join(("./configure", option)),
            shell=True,cwd=htslibLocalPath)
        if retcode != 0:
            return False
        else:
            return True
    except OSError as e:
        return False

# run make
def runMake():
    sys.stdout.flush()

    try:
        retcode = subprocess.call(
            ["make"],
            cwd=htslibLocalPath)
        if retcode != 0:
            return False
        else:
            return True
    except OSError as e:
        return False

# download and build htslib
def buildHtslib():
    # is downloaded?
    if not os.path.exists(htslibLocalPath):
        if not getHtslib():
            raise Exception("Could not download htslib")
        if not unpackHtslib():
            raise Exception("Could not unpack htslib")
    # is built?
    elif checkForHtslibSharedLibraries(htslibLocalPath):
        if not runConfigure(option=""):
            raise Exception("Could not run htslib configure")
        if not runMake():
            raise Exception("Could not run htslib make")

# check for htslib shared library files
def checkForHtslibSharedLibraries(dir):
    # check for libhts.dylib or libhts.so
    if not os.path.exists(os.path.join(dir,"libhts" + SHARED_EXT)):
        print(os.path.join(dir,"libhts" + SHARED_EXT))
        return False
    # check for libhts.3.dylib
    if not os.path.exists(os.path.join(dir,"libhts" + ".3" + SHARED_EXT)):
        # check for libhts.so.3
        if not os.path.exists(os.path.join(dir,"libhts" + SHARED_EXT + ".3")):
            sys.stderr.write("Installed htslib is too out of date (needed version >= 1.10)\n")
            return False
    return True
        
# get HTSLIB_DIR env var
HTSLIB_DIR = os.environ.get("HTSLIB_DIR",None)

# Find htslib under /usr/local/lib
# or specify it with HTSLIB_DIR env var
# or download and build it
def resolveHtslib():
    if HTSLIB_DIR:
        if not checkForHtslibSharedLibraries(HTSLIB_DIR):
            raise Exception("No htslib usable shared library in HTSLIB_DIR")
        else:
            return os.path.abspath(HTSLIB_DIR)
    elif checkForHtslibSharedLibraries("/usr/local/lib/"):
        return "/usr/local/lib/"
    else:
        sys.stderr.write("Downloading and building htslib {}\n".format(htslibVersion))
        buildHtslib()
        if not checkForHtslibSharedLibraries(htslibLocalPath):
            raise Exception("Something went wrong during the build process: no shared library files")
        return htslibLocalPath

class MyCleaner(clean):

    # cleans up htslib install 
    # then runs normal clean 
    def run(self):
        if os.path.exists(htslibLocalPath):
            shutil.rmtree(htslibLocalPath)
        if os.path.exists("htslib-{v}.tar.bz2".format(v=htslibVersion)):
            os.remove("htslib-{v}.tar.bz2".format(v=htslibVersion))
        super().run()

htslib_shared_path = resolveHtslib()
recontigSources = glob.glob(os.path.join("source","recontig","*.d"))
dhtslibSources = glob.glob(os.path.join("dhtslib","source","dhtslib","*.d")) + \
                glob.glob(os.path.join("dhtslib","source","dhtslib", "bed","*.d")) +  \
                glob.glob(os.path.join("dhtslib","source","dhtslib", "gff","*.d")) +  \
                glob.glob(os.path.join("dhtslib","source","dhtslib", "sam","*.d")) +  \
                glob.glob(os.path.join("dhtslib","source","dhtslib", "vcf","*.d")) +  \
                glob.glob(os.path.join("dhtslib","source","htslib","*.d"))

setup(
    name=projName,
    version='1.0.0',
    platforms=["POSIX", "UNIX", "MacOS"],
    # packages=['htslib-{}'.format(htslibVersion)],
    # package_data={'htslib-{}'.format(htslibVersion): ['htslib-{}/libhts.so'.format(htslibVersion)]},
    ext_modules=[
        Extension(projName, recontigSources + dhtslibSources,
            extra_compile_args=['-w','-L-lhts'],
            build_deimos=True,
            d_lump=True,
            library_dirs=[htslib_shared_path],
            libraries = ["hts"]
        ),
    ],
    cmdclass={"clean": MyCleaner}
)