########用于从原始二级谱图中识别指定肽段的b,y离子
import numpy as np
import pandas as pd
from constant import Mass

class mz_cal():   ####计算肽和交联肽一级二级m/z的类
    def extract_from_cl_peptide(self, cl_peptide):  ####用于从plink的交联肽序列中提取两条肽和位点的信息
        id1 = cl_peptide.find('X')
        pep1 = cl_peptide[:id1]
        pep2 = cl_peptide[(id1 + 1):]
        site1 = pep1.find('U')
        site2 = pep2.find('U')
        peptide1 = pep1.replace('U', '')
        peptide2 = pep2.replace('U', '')
        return peptide1, peptide2, site1, site2

    def regular_peptide_m_z(self, peptide, charge):
        mass = 0
        for i in range(len(peptide)):
            mass = mass + Mass.AA_residue_mass[peptide[i]]
        mass = mass + Mass.loss_mass['H2O']
        m_z = mass / charge + Mass.loss_mass['H+']
        return m_z

    def crosslink_peptide_m_z(self, cl_peptide, charge, crosslinker):
        peptide1, peptide2, site1, site2 = self.extract_from_cl_peptide(cl_peptide)
        mass1 = self.calculate_regular_peptide_m_z(peptide1, 1) - Mass.loss_mass['H+']
        mass2 = self.calculate_regular_peptide_m_z(peptide2, 1) - Mass.loss_mass['H+']
        mass = Mass.crosslinker_mass[crosslinker] + mass1 + mass2
        m_z = mass / charge + Mass.loss_mass['H+']
        return m_z

    def regular_peptide_msms_m_z(self, peptide, by_type, num, charge, loss):  #####计算普通肽的二级b，y离子的m/z，
        if charge > 2:  ###对普通肽而言，大于2价的碎片不考虑
            return -1
        else:
            if by_type == 'b':
                b = peptide
            elif by_type == 'y':
                b = peptide[::-1]
            mass = 0
            for j in range(num):
                mass = mass + Mass.AA_residue_mass[b[j]]
            mass = mass - Mass.loss_mass[loss]
            if by_type == 'b':
                mass = mass + 1.00783 ###加上氢原子
                m_z = mass / charge + (charge - 1)*1.00728/charge  ###因为碎裂后C原子带正点，所以对于带N个正电的时候，只需要加N-1个H+
            elif by_type == 'y':
                mass = mass + 17.00274  ###加上羟基
                m_z = mass / charge + (charge + 1)*1.00728/charge   ###碎裂后带有原来肽的H+，但是和N上的负电中和了，所以要带N个正点的话，需要加N+1个H+
            if (m_z < 200) or (m_z > 1700):
                return -1
            else:
                return m_z

    def crosslink_peptide_msms_m_z(self, cl_peptide, crosslinker, by_type, num, charge, loss):
        peptide1, peptide2, site1, site2 = self.extract_from_cl_peptide(cl_peptide)
        peptide1_mass = 0
        peptide2_mass = 0
        mass = 0
        m = 1
        for i in range(len(peptide1)):
            peptide1_mass = peptide1_mass + Mass.AA_residue_mass[peptide1[i]]
        for i in range(len(peptide2)):
            peptide2_mass = peptide2_mass + Mass.AA_residue_mass[peptide2[i]]
        peptide1_mass = peptide1_mass + Mass.loss_mass['H2O']
        peptide2_mass = peptide2_mass + Mass.loss_mass['H2O']
        if by_type == '1b':
            a = peptide1
            if num >= len(a):
                m = -1
            else:
                for i in range(num):
                    mass = mass + Mass.AA_residue_mass[a[i]]
                if num >= site1:
                    mass = mass + peptide2_mass + Mass.crosslinker_mass[crosslinker]
                mass = mass + 1.00783  ###加上氢原子
                mass = mass - Mass.loss_mass[loss]
                m_z = mass / charge + (charge - 1) * 1.00728 / charge  ###因为碎裂后C原子带正点，所以对于带N个正电的时候，只需要加N-1个H+
        elif by_type == '1y':
            a = peptide1[::-1]
            if num >= len(a):
                m = -1
            else:
                for i in range(num):
                    mass = mass + Mass.AA_residue_mass[a[i]]
                if num >= (len(peptide1) - site1 + 1):
                    mass = mass + peptide2_mass + Mass.crosslinker_mass[crosslinker]
                mass = mass + 17.00274  ###加上羟基
                mass = mass - Mass.loss_mass[loss]
                m_z = mass / charge + (charge + 1) * 1.00728 / charge  ###碎裂后带有原来肽的H+，但是和N上的负电中和了，所以要带N个正点的话，需要加N+1个H+
        elif by_type == '2b':
            a = peptide2
            if num >= len(a):
                m = -1
            else:
                for i in range(num):
                    mass = mass + Mass.AA_residue_mass[a[i]]
                if num >= site2:
                    mass = mass + peptide1_mass + Mass.crosslinker_mass[crosslinker]
                mass = mass + 1.00783  ###加上氢原子
                mass = mass - Mass.loss_mass[loss]
                m_z = mass / charge + (charge - 1) * 1.00728 / charge  ###因为碎裂后C原子带正点，所以对于带N个正电的时候，只需要加N-1个H+
        elif by_type == '2y':
            a = peptide2[::-1]
            if num >= len(a):
                m = -1
            else:
                for i in range(num):
                    mass = mass + Mass.AA_residue_mass[a[i]]
                if num >= (len(peptide2) - site2 + 1):
                    mass = mass + peptide1_mass + Mass.crosslinker_mass[crosslinker]
                mass = mass + 17.00274  ###加上羟基
                mass = mass - Mass.loss_mass[loss]
                m_z = mass / charge + (charge + 1) * 1.00728 / charge  ###碎裂后带有原来肽的H+，但是和N上的负电中和了，所以要带N个正点的话，需要加N+1个H+
        if m == -1:
            m_z = -1
        if (m_z < 200) or (m_z > 1700):
            m_z = -1
        return m_z



