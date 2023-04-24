import numpy as np
import pandas as pd
import argparse

def get_args():  ##设置需要传入的参数
    parser = argparse.ArgumentParser(description='convert bruker mgf to plink')
    parser.add_argument('--filedir', type=str, default=None)  ###输入文件的路径
    return parser.parse_args()

def filter_peptide_with_score_in_msms(datadir):   ####根据电荷态和PSM的Qvalue来filter
    data = pd.read_csv(datadir)
    data1 = data.sort_values('combine_peptide_z', ignore_index=True)  ####按照肽名字排序
    peptide_list = np.array(data1['combine_peptide_z'])
    name = list(data1)
    data5 = np.array(data1)  ###用已经有的矩阵来存放新产生的数据，这样可以大大提高速度
    peptide = peptide_list[0]
    index_list = []
    peptide_num = len(set(data1['combine_peptide_z']))
    num = 0
    lenth1 = 0
    for i in range(len(peptide_list)):
        if peptide_list[i] == peptide:
            index_list.append(i)
        else:
            num = num + 1
            print(num,'in',peptide_num)
            data2 = data1.iloc[index_list]
            q_list = list(data2['score'])
            file_list = list(data2['filename'])
            psm_list = list(data2['title'])
            file = file_list[int(q_list.index(min(q_list)))]
            psm = psm_list[int(q_list.index(min(q_list)))]
            data3 = data2[data2['filename'] == file]
            # data4 = data3[data3['title'] == psm]
            lenth2 = len(data3['combine_peptide_z'])
            data5[lenth1:(lenth1+lenth2),:] = np.array(data3)  ###把filter出来的数据存入data5
            lenth1 = lenth1 + lenth2
            index_list = []
            index_list.append(i)
            peptide = peptide_list[i]
    data2 = data1.iloc[index_list]
    q_list = list(data2['score'])
    file_list = list(data2['filename'])
    psm_list = list(data2['title'])
    file = file_list[int(q_list.index(min(q_list)))]
    psm = psm_list[int(q_list.index(min(q_list)))]
    data3 = data2[data2['filename'] == file]
    # data4 = data3[data3['title'] == psm]
    lenth2 = len(data3['combine_peptide_z'])
    data5[lenth1:(lenth1 + lenth2), :] = np.array(data3)
    lenth1 = lenth1 + lenth2
    data6 = pd.DataFrame(data5[:lenth1,:],columns=name)
    outputdir = datadir.split('.csv')[0] + '_filter.csv'
    data6.to_csv(outputdir,index = False)
    data7 = pd.DataFrame()
    data7['ModifiedPeptide'] = data6['combine_peptide']
    data7['PrecursorCharge'] = data6['charge']
    data7['PrecursorMz'] = data6['m_z']
    data7['FragmentCharge'] = data6['Fragment_charge']
    data7['ProductMz'] = data6['Fragment_m_z_calculation']
    data7['Tr_recalibrated'] = data6['rt']
    data7['IonMobility'] = data6['k0']
    data7['LibraryIntensity'] = data6['Fragment_intensity']
    outputdir1 = datadir.split('.csv')[0] + '_filter_DIANN_lib.csv'
    data7.to_csv(outputdir1, index=False)

if __name__ == '__main__':
    args = get_args()
    filter_peptide_with_score_in_msms(args.filedir)