import pandas as pd
import numpy as np
import math
import time
import argparse
import shutil
import os
######本程序用于将plink生成的结果转换成常规的csv表格并提取其中的rt，k0等关键信息，方便统计分析
def get_args():  ##设置需要传入的参数
    parser = argparse.ArgumentParser(description='Predict msms of peptides')
    parser.add_argument('--inputdir', type=str, default=None)  ###plink的report文件夹中交联肽peptide的csv文件名
    parser.add_argument('--outputdir', type=str, default='file')  ###输出的文件路径不包含后缀
    return parser.parse_args()


def extract_from_symbol(str, symbol):  ###定义一个函数来提取相同的开始符和结束符中间的字符串,symbol可以设置为任意字符
    index = []
    for i in range(len(str)):
        if str[i] == symbol:
            index.append(i)
    start_index = index[0]
    end_index = index[1]
    str1 = str[(start_index + 1):end_index]
    return str1

def extract_from_cl_peptide(cl_peptide):  ####用于从plink的交联肽序列中提取两条肽和位点的信息
    id1 = cl_peptide.find('-')
    pep1 = cl_peptide[:id1]
    pep2 = cl_peptide[(id1+1):]
    id1 = pep1.find('(')
    id2 = pep1.find(')')
    peptide1 = pep1[:id1]
    site1 = pep1[(id1+1):id2]
    id1 = pep2.find('(')
    id2 = pep2.find(')')
    peptide2 = pep2[:id1]
    site2 = pep2[(id1 + 1):id2]
    return peptide1, peptide2, site1, site2

def extract_from_mono_peptide(mono_peptide):  ####用于从plink的交联肽序列中提取两条肽和位点的信息
    id1 = mono_peptide.find('(')
    id2 = mono_peptide.find(')')
    peptide = mono_peptide[:id1]
    site = mono_peptide[(id1+1):id2]
    return peptide, site

def extract_from_loop_peptide(loop_peptide):  ####用于从plink的交联肽序列中提取两条肽和位点的信息
    lenth = len(loop_peptide)
    id1 = loop_peptide.find('(')
    id2 = loop_peptide.find(')')
    peptide = loop_peptide[:id1]
    site1 = loop_peptide[(id1+1):id2]
    site2 = loop_peptide[(id2+2):(lenth-1)]
    return peptide, site1, site2

def calculate_ccs(peptide_m_z, peptide_charge, peptide_k0):
    m = 28.00615
    t = 304.7527
    coeff = 18500 * peptide_charge * math.sqrt((peptide_m_z * peptide_charge + m) / (peptide_m_z * peptide_charge * m * t))
    ccs = coeff*peptide_k0
    return ccs

class mgf_from_bruker():  ######适用于brukerDA软件生成的mgf文件经过我们自己的程序转换后产生的mgf文件的plink2搜库结果
    def change_plink_filter_crosslink(self, plinkfile_dir, results_dir):  ####这里将plink生成的filtered_cross-linked_peptides文件转换为常规的文件形式方便统计, plinkfile_dir是要读取的plink文件路径，results_dir是结果文件保存路径
        plinkfile_dir1 = plinkfile_dir
        plinkfile_dir2 = plinkfile_dir.split('.csv')[0] + '_tmp.csv'  ###生成一个复制的临时文件用于添加第一行，方便使用pandas读取
        shutil.copy(plinkfile_dir1, plinkfile_dir2)
        with open(plinkfile_dir2, "r+") as f:
            contents = f.read()
            f.seek(0, 0)
            f.write('1,2,3,4,5,6,7,8,9,10,11' + '\n' + contents)  ###在文件第一行写入十一列的行名，方便pandas读取CSV文件
        data = pd.read_csv(plinkfile_dir2)
        data = np.array(data)
        cl_peptide_list, m_z_list, charge_list, peptide1_list, site1_list, peptide2_list, site2_list, rt_list, k0_list, ccs_list, score_list, precursor_Mass_Error_list, intensity_list, type_list = [], [], [], [], [], [], [], [], [], [], [], [], [], []
        title_list = []
        filename_list = []
        cmpd_list = []
        for i in range(2,len(data)):
            if str(data[i,0]) != 'nan':
                peptide = data[i,1]
                peptide1, peptide2, site1, site2 = extract_from_cl_peptide(peptide)
                cl_peptide = data[i,1]
                modif = str(data[i,3])
                peptide1, peptide2 = modif_clpeptide(peptide1, peptide2, modif)
                mass = data[i,2]
                type = 3
            if str(data[i,0]) == 'nan':
                if ('U' in str(peptide)) or ('O' in str(peptide)) or ('X' in str(peptide)) or ('B' in str(peptide)) or ('J' in str(peptide)) or ('Z' in str(peptide)):  ###删除带有U,O的肽
                    print(1)
                else:
                    title = data[i,2]
                    a = title.index('.')
                    filename = title[:a]
                    rt = extract_from_symbol(title, '$')
                    k0 = extract_from_symbol(title, '#')
                    intensity = extract_from_symbol(title, '|')
                    find_all = lambda c, s: [x for x in range(c.find(s), len(c)) if c[x] == s]  ###定义函数，寻找某个字符s的所有index
                    id1 = find_all(title, '.')
                    id2 = find_all(title, '.')
                    cmpd = title[(id1[5] + 1):id2[6]]
                    charge = data[i,3]
                    score = data[i,6]
                    precursor_Mass_Error = data[i,8]
                    m_z = m_z = (float(mass) + (int(charge) - 1) * 1.00728)/ int(charge)
                    ccs = calculate_ccs(m_z, int(charge), float(k0))
                    cl_peptide_list.append(cl_peptide)
                    m_z_list.append(m_z)
                    charge_list.append(charge)
                    peptide1_list.append(peptide1)
                    peptide2_list.append(peptide2)
                    site1_list.append(site1)
                    site2_list.append(site2)
                    rt_list.append(rt)
                    k0_list.append(k0)
                    ccs_list.append(ccs)
                    score_list.append(score)
                    precursor_Mass_Error_list.append(precursor_Mass_Error)
                    title_list.append(title)
                    filename_list.append(filename)
                    intensity_list.append(intensity)
                    cmpd_list.append(cmpd)
                    type_list.append(type)
        data = pd.DataFrame()
        data['title'] = title_list
        data['filename'] = filename_list
        data['peptide'] = cl_peptide_list
        data['m_z'] = m_z_list
        data['charge'] = charge_list
        data['peptide1'] = peptide1_list
        data['peptide2'] = peptide2_list
        data['site1'] = site1_list
        data['site2'] = site2_list
        data['rt'] = rt_list
        data['k0'] = k0_list
        data['ccs'] = ccs_list
        data['score'] = score_list
        data['precursor_Mass_Error(ppm)'] = precursor_Mass_Error_list
        data['intensity'] = intensity_list
        data['cmpd'] = cmpd_list
        data['peptide_type'] = type_list
        combine_pep, combine_pep_z = [], []
        for i in range(len(peptide1_list)):
            pep1 = peptide1_list[i]
            pep2 = peptide2_list[i]
            s1 = int(site1_list[i])
            s2 = int(site2_list[i])
            z = int(charge_list[i])
            pep1 = list(pep1)
            pep2 = list(pep2)
            pep1.insert(s1, 'U')
            pep2.insert(s2, 'U')
            pep1 = ''.join(pep1)
            pep2 = ''.join(pep2)
            pep = pep1 + 'X' + pep2
            pep_z = pep1 + 'X' + pep2 + str(z)
            combine_pep.append(pep)
            combine_pep_z.append(pep_z)
        data['combine_peptide'] = combine_pep
        data['combine_peptide_z'] = combine_pep_z
        os.remove(plinkfile_dir2)
        results_dir = f'{results_dir}_crosslink.csv'
        data.to_csv(results_dir, index=False)


def modif_clpeptide(pep1,pep2,mod):
    len1 = len(pep1)
    len2 = len(pep2)
    if 'M' in mod:
        find_all = lambda c, s: [x for x in range(c.find(s), len(c)) if c[x] == s]  ###定义函数，寻找某个字符s的所有index
        id_list = find_all(mod, 'M')
        for id in id_list:
            if mod[(id + 4)] == ')':
                modid = int(mod[(id + 3)])
                if modid < (len1 + len2 + 1):
                    if modid > len1:
                        b = list(pep2)
                        if b[(modid - len1 - 1)] == 'M':
                            b[(modid - len1 - 1)] = 'B'
                        pep2 = ''.join(b)
                    else:
                        b = list(pep1)
                        if b[(modid - len1 - 1)] == 'M':
                            b[(modid - len1 - 1)] = 'B'
                        pep1 = ''.join(b)
            elif mod[(id + 5)] == ')':
                modid = int(mod[(id + 3):(id + 5)])
                if modid < (len1 + len2 + 1):
                    if modid > len1:
                        b = list(pep2)
                        if b[(modid - len1 - 1)] == 'M':
                            b[(modid - len1 - 1)] = 'B'
                        pep2 = ''.join(b)
                    else:
                        b = list(pep1)
                        if b[(modid - len1 - 1)] == 'M':
                            b[(modid - len1 - 1)] = 'B'
                        pep1 = ''.join(b)
    return pep1, pep2

def modif_singlepeptide(pep,mod):
    if 'M' in mod:
        find_all = lambda c, s: [x for x in range(c.find(s), len(c)) if c[x] == s]  ###定义函数，寻找某个字符s的所有index
        id_list = find_all(mod, 'M')
        for id in id_list:
            if mod[(id + 4)] == ')':
                modid = int(mod[(id + 3)])
                b = list(pep)
                if b[(modid - 1)] == 'M':
                    b[(modid - 1)] = 'B'
                pep = ''.join(b)
            elif mod[(id + 5)] == ')':
                modid = int(mod[(id + 3):(id + 5)])
                b = list(pep)
                if b[(modid - 1)] == 'M':
                    b[(modid - 1)] = 'B'
                pep = ''.join(b)
    return pep

def extract_highest_intensity(inputdir, outputdir):
    data = pd.read_csv(inputdir)
    data_filter_index = 0
    name = list(data)
    for i in [3,4,5]:
        data1 = data[data['charge']==i]
        peptide = list(data1['peptide'])
        uni_pep = set(peptide)
        for pep in uni_pep:
            data2 = data1[data1['peptide']==pep]
            intensity = list(data2['intensity'])
            data3 = data2[data2['intensity']==max(intensity)]
            data3 = np.array(data3)
            if data_filter_index == 0:
                data_filter = data3[0,]
                data_filter_index = 1
            else:
                data_filter = np.row_stack((data_filter,data3[0,]))
    data_filter = pd.DataFrame(data_filter,columns=name)
    data_filter.to_csv(outputdir,index=False)

def extract_lowest_score(inputdir, outputdir):
    data = pd.read_csv(inputdir)
    data_filter_index = 0
    name = list(data)
    for i in [3,4,5]:
        data1 = data[data['charge']==i]
        peptide = list(data1['peptide'])
        uni_pep = set(peptide)
        for pep in uni_pep:
            data2 = data1[data1['peptide']==pep]
            intensity = list(data2['score'])
            data3 = data2[data2['score']==min(intensity)]
            data3 = np.array(data3)
            if data_filter_index == 0:
                data_filter = data3[0,]
                data_filter_index = 1
            else:
                data_filter = np.row_stack((data_filter,data3[0,]))
    data_filter = pd.DataFrame(data_filter,columns=name)
    data_filter.to_csv(outputdir,index=False)

def extract_from_peptide(filedir):
    data = pd.read_csv(filedir)
    peptide = data['peptide']
    peptide1, peptide2, site1, site2, type_list = [], [], [], [], []
    for pep in peptide:
        if '-' in pep:
            pep1, pep2, s1, s2 = extract_from_cl_peptide(pep)
            type = 3
        else:
            if '(' in pep:
                b = list(pep)
                num = b.count('(')
                if num == 1:
                    pep1, s1 = extract_from_mono_peptide(pep)
                    pep2, s2 = -1, -1
                    type = 1
                elif num == 2:
                    pep1, s1, s2 = extract_from_loop_peptide(pep)
                    pep2 = -1
                    type = 2
            else:
                pep1 = pep
                pep2, s1, s2 = -1, -1, -1
                type = 0
        peptide1.append(pep1)
        peptide2.append(pep2)
        site1.append(s1)
        site2.append(s2)
        type_list.append(type)
    data['peptide1'] = peptide1
    data['peptide2'] = peptide2
    data['site1'] = site1
    data['site2'] = site2
    data['type'] = type_list
    data.to_csv(filedir,index=False)

def filter_plink_results(datadir):   ####筛选相同肽中score最小的peptide
    data = pd.read_csv(f'{datadir}.csv')
    charge_list = [3, 4, 5]
    x = 0
    for z in charge_list:
        data1 = data[data['charge'] == z]
        data1 = data1.sort_values('peptide', ignore_index=True)  ####按照肽名字排序
        peptide_list = np.array(data1['peptide'])
        name = list(data1)
        data5 = np.array(data1)  ###用已经有的矩阵来存放新产生的数据，这样可以大大提高速度
        peptide = peptide_list[0]
        index_list = []
        peptide_num = len(set(data1['peptide']))
        num = 0
        lenth1 = 0
        for i in range(len(peptide_list)):
            if peptide_list[i] == peptide:
                index_list.append(i)
            else:
                num = num + 1
                print(num,'in',peptide_num)
                aa = time.clock()
                data2 = data1.iloc[index_list]
                bb = time.clock()
                print(bb-aa)
                aa = time.clock()
                q_list = list(data2['score'])
                psm_list = list(data2['title'])
                psm = psm_list[int(q_list.index(min(q_list)))]
                data3 = data2[data2['title'] == psm]
                lenth2 = len(data3['peptide'])
                data5[lenth1:(lenth1+lenth2),:] = np.array(data3)  ###把filter出来的数据存入data5
                lenth1 = lenth1 + lenth2
                index_list = []
                index_list.append(i)
                peptide = peptide_list[i]
                bb = time.clock()
                print(bb - aa)
        data2 = data1.iloc[index_list]
        q_list = list(data2['score'])
        psm_list = list(data2['title'])
        psm = psm_list[int(q_list.index(min(q_list)))]
        data3 = data2[data2['title'] == psm]
        lenth2 = len(data3['peptide'])
        data5[lenth1:(lenth1 + lenth2), :] = np.array(data3)  ###把filter出来的数据存入data5
        lenth1 = lenth1 + lenth2
        data6 = pd.DataFrame(data5[:lenth1,:],columns=name)
        if x == 0:
            data7 = data6
            x = 1
        else:
            data7 = np.row_stack((data7, data6))
    data7 = pd.DataFrame(data7, columns=name)
    data7.to_csv(f'{datadir}_charge345_ion.csv',index = False)

def process_mgf_from_bruker(plinkfile, resultdir):
    x = mgf_from_bruker()
    x.change_plink_filter_crosslink(plinkfile, resultdir)
    extract_lowest_score(f'{resultdir}_crosslink.csv', f'{resultdir}_crosslink_filter.csv')

if __name__ == '__main__':
    args = get_args()
    process_mgf_from_bruker(args.inputdir, args.outputdir)





