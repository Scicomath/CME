# 重离重离子碰撞中背景磁场的计算
本课题是我的本科毕业论文，主要讨论重离重离子碰撞中背景磁场的计算以及背景磁场的性质。

## 代码运行说明
目前代码都是用C语言写的，其中的积分程序用到了GSL(GNU Scientific Library)库。所以要先安装GSL才能运行。安装步骤如下：
- 到[GSL官网](http://www.gnu.org/software/gsl/)下载最新版的安装包
- 将下载的压缩包解压：`$ tar -xvf gsl-latest.tar.gz`
- 切换到解压后的目录中：`$ cd gsl-1.16`
- 确保自己安装了编译器，然后运行：`$ ./configure`
- `$ make`
- `$ make install`
- 下面一步很关键，现在直接运行一般会提示找不到XXX.so。这是因为系统没有到正确到位置找它。所以要编辑/etc/ld.so.conf文件，在后面加入"/usr/local/lib"这一行，保存之后，再运行：`$ /sbin/ldconfig -v`更新一下配置就可以啦。