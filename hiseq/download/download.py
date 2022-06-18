#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import os
import sys
import yaml 
# import oss2
from hiseq.utils.file import file_abspath, file_exists, check_path


# parse config from txt file
def parse_oss(x, out=None):
    """
    parse the email content for AccessKeyId, AccessKeySecret, 预设OSS路径
    also for:
    合同号, 开题单号
    """
    # CN to EN
    k = {
        '合同号': 'contract_id',
        '开题单号': 'project_id',
        'AccessKeyId': 'AccessKeyId',
        'AccessKeySecret': 'AccessKeySecret',
        '预设OSS路径': 'OSS_path',
        '过期时间': 'expire',
    }
    d = {'endpoint': 'http://oss-cn-shanghai.aliyuncs.com'}
    try:
        with open(x) as r:
            for l in r:
                s = l.strip().split(None, 1)
                if len(s) == 2:
                    if s[0] in k:
                        kk = k.get(s[0])
                        d.update({kk:s[1]})
    except:
        print('failed reading file: {}'.format(x))
        d = None
    # save as yaml
    try:
        with open(out, 'wt') as w:
            yaml.dump(d, w)
    except:
        print('write to file, failed: {}'.format(out))
    return d


# use ossutil commandline tool
def download(x, outdir):
    """
    Download OSS files using ossutils 
    """
    # 1. load config files
    x = file_abspath(x)
    c = os.path.join(os.path.dirname(x), 'config.yaml')
    d = parse_oss(x, c)
    print(x, c, d)
    # 2. prep command
    outdir = file_abspath(outdir)
    check_path(outdir)
    ossutil = '/data/biosoft/ossutil_aliyuncs/ossutil64'
    cmd = ' '.join([
        ossutil,
        'sync',
        '-e {}'.format(d.get('endpoint')),
        '-i {}'.format(d.get('AccessKeyId')),
        '-k {}'.format(d.get('AccessKeySecret')),
        d.get('OSS_path'),
        outdir
    ])
    cmd_txt = os.path.join(outdir, 'run.sh')
    with open(cmd_txt, 'wt') as w:
        w.write(cmd+'\n')
    # 3. run
    os.system(cmd)


def main():
    if len(sys.argv) < 3:
        print('Usage: python download.py <oss.txt> <outdir>')
        sys.exit(1)
    download(sys.argv[1], sys.argv[2])
    
    
if __name__ == '__main__':
    main()
    
# EOF


# download file/folder from oss using SDK (oss2) 
# SDK：提供丰富、完整的各类语言SDK demo，易于开发。
# 上传文件夹：SDK不支持直接上传文件夹，您可以在上传时设置相同的文件名前缀，
# 并使用正斜线（/）隔开。例如上传a.txt、b.txt、c.txt三个文件到abc文件夹，
# 在上传时设置ObjectName为abc/a.txt、abc/b.txt、abc/c.txt即可。
# 下载文件夹：SDK不支持直接下载文件夹，您可以在下载时将文件下载到同一个本地文件夹中。

# 断点下载。
# 实现的方法是：
# 在本地创建一个临时文件，文件名由原始文件名加上一个随机的后缀组成；
# 通过指定请求的 Range 头按照范围并发读取OSS文件，并写入到临时文件里对应的位置；
# 全部完成之后，把临时文件重命名为目标文件 （即 filename ）
# 在上述过程中，断点信息，即已经完成的范围，会保存在磁盘上。因为某种原因下载中断，后续如果下载 同样的文件，也就是源文件和目标文件一样，就会先读取断点信息，然后只下载缺失的部分。
# 缺省设置下，断点信息保存在 HOME 目录的一个子目录下。可以通过 store 参数更改保存位置。
# 使用该函数应注意如下细节：
# 对同样的源文件、目标文件，避免多个程序（线程）同时调用该函数。因为断点信息会在磁盘上互相覆盖，或临时文件名会冲突。
# 避免使用太小的范围（分片），即 part_size 不宜过小，建议大于或等于 oss2.defaults.multiget_part_size 。
# 如果目标文件已经存在，那么该函数会覆盖此文件。

# 以下代码展示了文件下载的用法，如下载文件、范围下载、断点续传下载等。
# 首先初始化AccessKeyId、AccessKeySecret、Endpoint等信息。
# 通过环境变量获取，或者把诸如“<你的AccessKeyId>”替换成真实的AccessKeyId等。
#
# 以杭州区域为例，Endpoint可以是：
#   http://oss-cn-hangzhou.aliyuncs.com
#   https://oss-cn-hangzhou.aliyuncs.com
# 分别以HTTP、HTTPS协议访问。
# access_key_id = os.getenv('OSS_TEST_ACCESS_KEY_ID', '<你的AccessKeyId>')
# access_key_secret = os.getenv('OSS_TEST_ACCESS_KEY_SECRET', '<你的AccessKeySecret>')
# bucket_name = os.getenv('OSS_TEST_BUCKET', '<你的Bucket>')
# endpoint = os.getenv('OSS_TEST_ENDPOINT', '<你的访问域名>')
# # 确认上面的参数都填写正确了
# for param in (access_key_id, access_key_secret, bucket_name, endpoint):
#     assert '<' not in param, '请设置参数：' + param


# # 创建Bucket对象，所有Object相关的接口都可以通过Bucket对象来进行
# bucket = oss2.Bucket(oss2.Auth(access_key_id, access_key_secret), endpoint, bucket_name)

# key = 'motto.txt'
# content = oss2.to_bytes('a' * 1024 * 1024)
# filename = 'download.txt'

# # 上传文件
# bucket.put_object(key, content, headers={'content-length': str(1024 * 1024)})

# """
# 文件下载
# """

# # 下载文件
# result = bucket.get_object(key)

# # 验证一下
# content_got = b''
# for chunk in result:
#     content_got += chunk
# assert content_got == content

# # 下载到本地文件
# result = bucket.get_object_to_file(key, filename)

# # 验证一下
# with open(filename, 'rb') as fileobj:
#     assert fileobj.read() == content

# """
# 范围下载
# """

# # 范围下载，如果指定的范围无效，则下载整个文件
# result = bucket.get_object(key, byte_range=(0, 1023))

# # 验证一下
# content_got = b''
# for chunk in result:
#     content_got += chunk
# assert content_got == oss2.to_bytes('a'*1024)


# # 范围下载到本地文件
# result = bucket.get_object_to_file(key, filename, byte_range=(1024, 2047))

# # 验证一下
# with open(filename, 'rb') as fileobj:
#     assert fileobj.read() == oss2.to_bytes('a'*1024)


# """
# 断点续传下载
# """

# # 断点续传下载
# oss2.resumable_download(bucket, key, filename,
#                         multiget_threshold=200*1024,
#                         part_size=100*1024,
#                         num_threads=3)

# # 验证一下
# with open(filename, 'rb') as fileobj:
#     assert fileobj.read() == content

# # 清理文件
# os.remove(filename)