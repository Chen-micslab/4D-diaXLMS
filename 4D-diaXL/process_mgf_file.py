import pandas as pd
import argparse
import sys
######本程序用于将bruker生成的mgf文件转换为plink可读取的mgf文件

def get_args():  ##设置需要传入的参数
    parser = argparse.ArgumentParser(description='convert bruker mgf to plink')
    parser.add_argument('--filedir', type=str, default=None)  ###输入文件的路径
    parser.add_argument('--filename', type=str, default='file')  ###给文件命名
    return parser.parse_args()

def bruker_to_plink_mgf(file_path, file_name):
    data = open(file_path, 'r')
    reaction = 0
    data1 = []
    total_num = 0
    for line in data:
        if line.startswith('BEGIN'):
            total_num = total_num + 1
            reaction = 1
            if total_num % 200 == 0:
                sys.stdout.write("Number of processed spectrums: {}   \r".format(total_num))
                sys.stdout.flush()
            continue
        elif line.startswith('END'):
            if reaction == 1:
                data1.append('END IONS')
                data1.append(' ')
                reaction = 0
                continue
        if line.startswith('TITLE'):  ###需要提取一级的cmpd，还有scan1和scan2
            find_all = lambda c, s: [x for x in range(c.find(s), len(c)) if c[x] == s]  ###定义函数，寻找某个字符s的所有index
            ###提取一级信息
            id1 = find_all(line, ' ')
            id2 = find_all(line, ',')
            cmpd = line[(id1[0]+1):id2[0]]  ###一级质谱的cmpd，也不知道啥意思
            ###提取二级信息
            id1 = find_all(line, '#')
            id2 = find_all(line, ',')
            id3 = find_all(line, '-')
            if id3 == []:
                scan1 = line[(id1[0] + 1):id2[-1]]
                scan2 = line[(id1[0] + 1):id2[-1]]
            else:
                scan1 = line[(id1[0] + 1):id3[0]]
                scan2 = line[(id1[0] + 1):id3[0]]
                # scan2 = line[(id3[0] + 1):id2[-1]]
            ###提取淌度和保留时间的信息
            id1 = find_all(line, '=')
            id2 = find_all(line, ',')
            k0 = line[(id1[1]+1):id2[5]]  ###提取淌度
            id1 = find_all(line, ',')
            rt = line[(id1[3]+2):(id1[4]-3)]  ###提取保留时间
        elif line.startswith('RTINSECONDS'):  ###无需提取
            pass
        elif line.startswith('RAWSCANS'):  ###无需提取
            pass
        elif line.startswith('PEPMASS'):  ###提取其中的m/z
            find_all = lambda c, s: [x for x in range(c.find(s), len(c)) if c[x] == s]  ###定义函数，寻找某个字符s的所有index
            id1 = find_all(line, '=')
            id2 = find_all(line, '	')  ####注意这里不是普通空格
            m_z = line[(id1[0]+1):id2[0]]
            len1 = len(line)
            intensity = line[(id2[0]+1):(len1-1)]
            reaction = 0
        elif line.startswith('CHARGE'):  ###提取charge
            find_all = lambda c, s: [x for x in range(c.find(s), len(c)) if c[x] == s]  ###定义函数，寻找某个字符s的所有index
            id1 = find_all(line, '=')
            id2 = find_all(line, '+')
            charge = line[(id1[0] + 1):id2[0]]
            title = 'TITLE='+f'{file_name}'+f'.{scan1}'+f'.{scan2}'+f'|{intensity}|'+f'${rt}$'+f'#{k0}#'+f'.{charge}'+f'.{cmpd}.dta' ###在title中加入淌度和保留时间
            data1.append('BEGIN IONS')
            data1.append(title)
            mass = f'PEPMASS={m_z}'
            data1.append(mass)
            z = f'CHARGE={charge}+'
            data1.append(z)
            reaction = 1
        else:
            if reaction == 1:
                find_all = lambda c, s: [x for x in range(c.find(s), len(c)) if c[x] == s]  ###定义函数，寻找某个字符s的所有index
                id = find_all(line, '	')
                msms = line[:id[-1]]
                data1.append(msms)
    data1 = pd.DataFrame(data1,columns=['xyz'])
    sys.stdout.write("Number of processed spectrums: {}   \r".format(total_num))
    sys.stdout.flush()
    return data1


if __name__ == '__main__':
    args = get_args()
    data = bruker_to_plink_mgf(args.filedir, args.filename)
    outputdir = args.filedir.split('.mgf')[0] + '_plink.mgf'
    data.to_csv(outputdir, header= False, index= False)