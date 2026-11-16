### 导入必要的包
# from dbm.ndbm import library
import re, os, sys, logging, time, datetime
from optparse import OptionParser, OptionGroup
import subprocess
import threading

### 检查python的版本，我们需要使用python3
if sys.version_info < (3, 0):
    print("Python Version error: please use phthon3")
    sys.exit(-1)

### 程序版本号
_version = 'v0.1'


# -----------------------------------------------------------------------------------
# --- S 参数设置模块
# -----------------------------------------------------------------------------------
def configOpt():
    """Init for option
    """
    usage = 'Usage: %prog [-f] [other option] [-h]'
    p = OptionParser(usage)
    ##basic options
    p.add_option('-i', '--infile', dest='infile', action='store', type='string', default = "sc.all.csv", help='csv file')
    p.add_option('-o', '--outfile', dest='outfile', action='store', type='string',default = "sc.loom", help='loom file')

    group = OptionGroup(p, "Preset options")
    ##preset options
    group.add_option('-O', '--outDir', dest='outDir', default='./', action='store', type='string', help='output directory', metavar="DIR")
    group.add_option('-L', '--logPrefix', dest='logPrefix', default='', action='store', type='string', help='log file prefix')
    group.add_option('-E', '--email', dest='email', default='none', action='store', type='string', help='email address, if you want get a email when this job is finished,default is no email', metavar="EMAIL")
    group.add_option('-Q', '--quiet', dest='quiet', default=False, action='store_true', help='do not print messages to stdout')
    group.add_option('-K', '--keepTemp', dest='keepTemp', default=False, action='store_true', help='keep temp dir')
    group.add_option('-T', '--test', dest='isTest', default=False, action='store_true', help='run this program for test')
    p.add_option_group(group)
    opt, args = p.parse_args()
    return (p, opt, args)


def listToString(x):
    """获得完整的命令
    """
    rVal = ''
    for a in x:
        rVal += a + ' '
    return rVal


## 解析参数
opt_parser, opt, args = configOpt()



os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv(opt.infile);
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create(opt.outfile,x.X.transpose(),row_attrs,col_attrs);
# 行为基因，列为细胞
# lp.create(opt.outfile,x.X,row_attrs,col_attrs);