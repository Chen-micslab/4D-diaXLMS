import numpy as np
import pandas as pd
from mass_cal import mz_cal as M
import argparse
import sys

def get_args():  ##设置需要传入的参数
    parser = argparse.ArgumentParser(description='Generate 4D crosslinking library')
    parser.add_argument('--resultsdir', type=str, default=None)  ###plink的report文件夹中交联肽peptide的csv文件名
    parser.add_argument('--mgfdir', type=str, default=None)  ###输出的文件路径不包含后缀
    parser.add_argument('--outputdir', type=str, default=None)  ###输出的文件路径不包含后缀
    parser.add_argument('--crosslinker', type=str, default='DSS')  ###输出的文件路径不包含后缀
    return parser.parse_args()

def match_msms(specturm, m_z, precursor, spe, id_s):  ###m_z是一个list，输出匹配到list中荷质比的intensity
    title = 'TITLE='+precursor
    id_x = np.where(spe == title)
    id = id_s[int(id_x[0])]
    msms = []
    intensity = []
    id = id + 3  #########如果是bruker转的mgf就是id+3
    for i in range(id,len(specturm)):
        if specturm[i] == 'END IONS':
            break
        else:
            a = str(specturm[i])
            id1 = a.index('	')
            # id1 = a.index(' ')
            msms.append(float(a[:id1]))
            intensity.append(float(a[(id1+1):]))
    msms = np.array(msms)
    intensity = np.array(intensity)
    m_z_1 = []
    inten = []
    for mz in m_z:
        if mz == -1:
            m_z_1.append(-1), inten.append(-1)
        else:
            msms1 = np.abs(msms - mz)/msms
            if np.min(msms1) <= 0.00002:  #### 设置最大质量偏移为20ppm
                m_z_1.append(msms[np.argmin(msms1)]), inten.append(intensity[np.argmin(msms1)])
            else:
                m_z_1.append(0), inten.append(0)
    return m_z_1, inten

def crosslink_ion_generation(peptide1, peptide2):  #包含1b,1y,2b,2y: +1,+2: -NH3,-H2O,noloss
    len1 = len(peptide1)
    len2 = len(peptide2)
    z = ['1', '2', '3', '4', '5']   ####碎片可能电荷
    l = ['noloss'] * 5   ####中性丢失
    by_1 = (['1b'] * 5 + ['1y'] * 5) * (len1 - 1)
    c = []  ###碎裂位点
    for i in range(len1 - 1):
        j = i + 1
        c = c + [j] * 10
    for i in range(len2 - 1):
        j = i + 1
        c = c + [j] * 10
    by_2 = (['2b'] * 5 + ['2y'] * 5) * (len2 - 1)
    by = by_1 + by_2
    z = z * 2 * (len1 + len2 - 2)
    l = l * 2 * (len1 + len2 - 2)
    c = np.array(c)    ###碎裂位点
    by = np.array(by)   ###by离子类型
    z = np.array(z)   ###碎片电荷数
    l = np.array(l)   ###中性丢失
    data = np.column_stack((c, by, z, l))
    return data

def genenrate_all_crosslink_fragment(datadir, mgf_dir, outputdir, crosslinker):
    sys.stdout.write("Loading file......\r")
    spectrum = np.array(pd.read_csv(mgf_dir,sep='!'))  ####为了让所有数据都读在一列中，选一个最不可能出现的分隔符
    spectrum = spectrum.flatten()
    spe = []
    id_s = []
    for i in range(len(spectrum)):
        if str(spectrum[i]).startswith('TIT'):
            spe.append(spectrum[i])
            id_s.append(i)
    spe = np.array(spe)
    id_s = np.array(id_s)
    data = pd.read_csv(datadir)
    name = list(data)
    pep1 = np.array(data['peptide1'])
    pep2 = np.array(data['peptide2'])
    num = 0
    for i in range(len(pep1)):
        num = num + len(pep1[i]) + len(pep2[i]) - 2
    num = num * 10
    data = np.array(data)
    a = np.array(list(data[1,:]) + [0,0,0,0,0,0,0], dtype=object)  ###加入  dtype=object，可以防止字符串被截断
    data1 = np.tile(a, [num, 1])
    num = 0
    name = name + ['Fragment_num', 'Fragment_type', 'Fragment_charge', 'Neutral_loss']
    sys.stdout.write("Generating library......\r")
    for i in range(len(pep1)):
        peptide1 = pep1[i]
        peptide2 = pep2[i]
        len1 = (len(pep1[i]) + len(pep2[i]) - 2) * 10
        b = data[i,:]
        data_y = np.tile(b, [len1, 1])
        data_x = crosslink_ion_generation(peptide1, peptide2)
        data_xy = np.column_stack((data_y, data_x))
        data_xy = pd.DataFrame(data_xy, columns=name)
        cl_peptide = data_xy['combine_peptide']
        by_type = data_xy['Fragment_type']
        Fragment_num = data_xy['Fragment_num']
        Fragment_charge = data_xy['Fragment_charge']
        Neutral_loss = data_xy['Neutral_loss']
        title = data_xy['title']
        m_z = []
        m = M()
        for j in range(len(cl_peptide)):
            mz = m.crosslink_peptide_msms_m_z(cl_peptide[j], crosslinker, by_type[j], int(Fragment_num[j]), int(Fragment_charge[j]), Neutral_loss[j])
            m_z.append(mz)
        m_z_1, inten = match_msms(spectrum, m_z, title[0], spe, id_s)
        if np.sum(np.array(inten)>0) > 2:
            inten = np.array(inten)/np.max(inten)
            m_z = np.array(m_z)
            m_z_1 = np.array(m_z_1)
            data_xy = np.array(data_xy)
            data_xy = np.column_stack((data_xy,  m_z, m_z_1, inten))
            charge = int(data_xy[0, 3])
            for k in range(len(data_xy)):
                if int(data_xy[k,-5]) > charge:
                    data_xy[k, -1], data_xy[k, -2], data_xy[k, -3] = -1, -1, -1
            data1[num:(num + len1), :] = data_xy
            num = num + len1
    name = name + ['Fragment_m_z_calculation', 'Fragment_m_z_experiment', 'Fragment_intensity']
    for l in range(len(data1)):
        if data1[l,-1] < 0:
            data1[l, -1] = -1
    data1 = pd.DataFrame(data1, columns = name)
    data1 = data1[data1['Fragment_m_z_experiment'] > 0]
    data1 = data1[data1['Fragment_m_z_experiment'] < 1600]
    data1 = data1[data1['Fragment_intensity'] > 0]
    data1.to_csv(f'{outputdir}_normal_lib.csv', index = False)
    data2 = pd.DataFrame()
    data2['ModifiedPeptide'] = data1['combine_peptide']
    data2['PrecursorCharge'] = data1['charge']
    data2['PrecursorMz'] = data1['m_z']
    data2['FragmentCharge'] = data1['Fragment_charge']
    data2['ProductMz'] = data1['Fragment_m_z_calculation']
    data2['Tr_recalibrated'] = data1['rt']
    data2['IonMobility'] = data1['k0']
    data2['LibraryIntensity'] = data1['Fragment_intensity']
    data2.to_csv(f'{outputdir}_DIANN_lib.csv', index = False)
    return data1


if __name__ == '__main__':
    args = get_args()
    if args.outputdir == None:
        outputdir = args.resultsdir.split('.csv')[0]
    else:
        outputdir = args.outputdir
    genenrate_all_crosslink_fragment(args.resultsdir, args.mgfdir, outputdir, args.crosslinker)

