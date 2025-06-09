#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基因型数据拆分工具 - 重新实现版本
基于对 split_breed_plink.exe 的逆向分析结果

主要功能：
1. 读取PLINK格式的基因型数据(.map/.ped文件)
2. 根据芯片位置对应表按品种拆分数据
3. 生成分离后的基因型文件和压缩包
4. 支持中芯一号V1PLUS芯片格式验证

作者：基于逆向工程重新实现
日期：2025年
"""

import os
import sys
import time
import zipfile
from collections import defaultdict
from datetime import datetime

class GenotypeDataSplitter:
    """基因型数据拆分器"""
    
    def __init__(self):
        self.error_status = False
        self.end_time = "2024-12-31"
        
        # 中芯一号V1PLUS芯片的标准SNP列表
        self.v1plus_snps = [
            "CNCB10000416", "CNCB10002887", "CNCB10004677", 
            "CNCB10006046", "CNCB10009510", "CNCB10009951", "CNCB10010848"
        ]
        
        # 数据存储
        self.exclude_set = set()
        self.map_data = {}  # {prefix: [(snp_id, chromosome, position, alleles), ...]}
        self.ped_data = {}  # {prefix: [(sample_id, family_id, father_id, mother_id, sex, phenotype, genotypes), ...]}
        self.split_table = {}  # {chip_id: (sample_name, breed)}
        
    def check_expiry(self):
        """检查软件有效期"""
        now = datetime.now()
        end_date = datetime.strptime(self.end_time, "%Y-%m-%d")
        
        if now > end_date:
            print("错误：该软件已过期，请联系作者更新软件！")
            print("请按任意键退出窗口")
            input()
            sys.exit(1)
    
    def load_exclude_list(self, filename="exclude_chipid.txt"):
        """加载排除的芯片ID列表"""
        if not os.path.exists(filename):
            return
        
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                for line in f:
                    chip_id = line.strip()
                    if chip_id:
                        self.exclude_set.add(chip_id)
            print(f"已加载 {len(self.exclude_set)} 个排除的芯片ID")
        except Exception as e:
            print(f"读取排除列表失败: {e}")
    
    def load_split_table(self, filename="split_plink.txt"):
        """加载芯片位置对应表"""
        if not os.path.exists(filename):
            print(f"错误：工作目录中不存在芯片位置对应表 {filename} 文件")
            self.error_status = True
            return
        
        if os.path.getsize(filename) == 0:
            print(f"错误：芯片位置对应表 {filename} 为空文件")
            self.error_status = True
            return
        
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    parts = line.strip().split('\t')
                    if len(parts) != 3:
                        print(f"错误：在芯片位置对应表中第 {line_num} 行含有 {len(parts)} 列内容，不是 3 列内容")
                        self.error_status = True
                        continue
                    
                    chip_id, sample_name, breed = parts
                    
                    # 检查重复芯片号
                    if chip_id in self.split_table:
                        print(f"错误：在芯片位置对应表中存在重复芯片号 {chip_id}")
                        self.error_status = True
                        continue
                    
                    self.split_table[chip_id] = (sample_name, breed)
            
            print(f"已加载 {len(self.split_table)} 个芯片位置对应关系")
            
        except Exception as e:
            print(f"读取芯片位置对应表失败: {e}")
            self.error_status = True
    
    def find_plink_files(self):
        """查找工作目录中的PLINK文件"""
        map_files = []
        ped_files = []
        
        for filename in os.listdir('.'):
            if filename.endswith('.map'):
                prefix = filename[:-4]
                ped_file = prefix + '.ped'
                if os.path.exists(ped_file):
                    map_files.append(filename)
                    ped_files.append(ped_file)
                else:
                    print(f"错误：工作目录中存在 {filename} ，却不存在相应的 {ped_file}")
                    self.error_status = True
            elif filename.endswith('.ped'):
                prefix = filename[:-4]
                map_file = prefix + '.map'
                if not os.path.exists(map_file):
                    print(f"错误：工作目录中存在 {filename} ，却不存在相应的 {map_file}")
                    self.error_status = True
        
        if not map_files:
            print("错误：工作目录中不存在 map 文件")
            self.error_status = True
            return []
        
        return list(zip(map_files, ped_files))
    
    def load_map_file(self, filename):
        """加载MAP文件"""
        map_data = []
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    parts = line.strip().split('\t')
                    if len(parts) != 4:
                        print(f"错误：{filename} 文件中第 {line_num} 行不是4列")
                        self.error_status = True
                        continue
                    
                    chromosome, snp_id, genetic_distance, position = parts
                    map_data.append((snp_id, chromosome, position, genetic_distance))
            
            return map_data
        except Exception as e:
            print(f"读取MAP文件 {filename} 失败: {e}")
            self.error_status = True
            return []
    
    def load_ped_file(self, filename, expected_snp_count):
        """加载PED文件"""
        ped_data = []
        expected_columns = 6 + 2 * expected_snp_count
        
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    parts = line.strip().split('\t')
                    if len(parts) != expected_columns:
                        print(f"错误：在 {filename} 文件中第 {line_num} 行的列数为 {len(parts)}， "
                              f"不符合要求（正常列数应该为 6+2×位点数）")
                        self.error_status = True
                        continue
                    
                    family_id = parts[0]
                    sample_id = parts[1]
                    father_id = parts[2]
                    mother_id = parts[3]
                    sex = parts[4]
                    phenotype = parts[5]
                    genotypes = parts[6:]
                    
                    ped_data.append((sample_id, family_id, father_id, mother_id, sex, phenotype, genotypes))
            
            return ped_data
        except Exception as e:
            print(f"读取PED文件 {filename} 失败: {e}")
            self.error_status = True
            return []

    def validate_map_consistency(self):
        """验证MAP文件的一致性"""
        if len(self.map_data) <= 1:
            return True

        # 检查所有MAP文件的SNP数量是否一致
        snp_counts = [len(data) for data in self.map_data.values()]
        if len(set(snp_counts)) > 1:
            print("错误：工作目录中疑似存在多款芯片，不同map文件行数不一致")
            self.error_status = True
            return False

        # 检查SNP顺序是否一致
        first_snps = None
        for prefix, data in self.map_data.items():
            current_snps = [snp[0] for snp in data]  # 提取SNP ID
            if first_snps is None:
                first_snps = current_snps
            elif first_snps != current_snps:
                print("错误：工作目录中多个map文件的SNP顺序不一致")
                self.error_status = True
                return False

        return True

    def validate_v1plus_chip(self):
        """验证是否为中芯一号V1PLUS芯片格式"""
        if not self.map_data:
            return False

        # 获取第一个MAP文件的SNP列表
        first_map = list(self.map_data.values())[0]
        snp_ids = [snp[0] for snp in first_map]

        # 检查是否包含V1PLUS芯片的标准SNP
        v1plus_found = all(snp in snp_ids for snp in self.v1plus_snps)

        if not v1plus_found:
            print("错误：map文件中的SNP顺序与中芯一号V1PLUS芯片下机的原始map文件不一致")
            self.error_status = True
            return False

        return True

    def validate_data_integrity(self):
        """验证数据完整性"""
        # 检查芯片位置对应表中的芯片号是否在PED数据中存在
        all_sample_ids = set()
        for ped_data in self.ped_data.values():
            for sample_data in ped_data:
                sample_id = sample_data[0]
                all_sample_ids.add(sample_id)

        # 检查重复的个体号
        sample_id_counts = defaultdict(int)
        for sample_id in all_sample_ids:
            sample_id_counts[sample_id] += 1

        for sample_id, count in sample_id_counts.items():
            if count > 1:
                print(f"错误：在所有plink文件中存在重复芯片号 {sample_id}")
                self.error_status = True

        # 检查芯片位置对应表中的芯片号是否都有基因型数据
        missing_genotypes = []

        for chip_id in self.split_table:
            if chip_id not in all_sample_ids:
                missing_genotypes.append(chip_id)

        if missing_genotypes:
            for chip_id in missing_genotypes:
                print(f"错误：在芯片位置对应表中存在没有基因型的芯片号 {chip_id}")
            self.error_status = True

        # 检查芯片位置对应表的行数与基因型个体总数是否一致
        if len(self.split_table) != len(all_sample_ids):
            print("错误：芯片位置对应表的行数与基因型个体总数不一致")
            self.error_status = True

        return not self.error_status

    def split_by_breed(self):
        """按品种拆分基因型数据"""
        if self.error_status:
            return

        # 按品种分组
        breed_groups = defaultdict(list)
        for chip_id, (sample_name, breed) in self.split_table.items():
            breed_groups[breed].append((chip_id, sample_name))

        # 获取时间戳
        now_time_str = time.strftime("%Y%m%d%H%M%S", time.localtime())

        # 为每个品种生成文件
        for breed, samples in breed_groups.items():
            if len(samples) < 100:
                print(f"警告：品种 {breed} 的基因型个体数目小于 100")

            # 生成输出文件名
            output_prefix = f"{breed}_{now_time_str}"

            self.generate_breed_files(breed, samples, output_prefix)

    def generate_breed_files(self, breed, samples, output_prefix):
        """为特定品种生成输出文件"""
        # 创建输出文件
        ped_filename = f"{output_prefix}.ped"
        map_filename = f"{output_prefix}.map"
        idmap_filename = f"{output_prefix}.idmap.txt"
        zip_filename = f"{output_prefix}.zip"

        try:
            # 获取第一个MAP文件作为模板
            template_map = list(self.map_data.values())[0]

            # 写入MAP文件
            with open(map_filename, 'w', encoding='utf-8') as map_file:
                for snp_data in template_map:
                    snp_id, chromosome, position, genetic_distance = snp_data
                    map_file.write(f"{chromosome}\t{snp_id}\t{genetic_distance}\t{position}\n")

            # 收集该品种的基因型数据
            breed_ped_data = []
            id_mapping = []

            for chip_id, sample_name in samples:
                # 在所有PED数据中查找该芯片ID
                found = False
                for prefix, ped_data in self.ped_data.items():
                    for sample_data in ped_data:
                        if sample_data[0] == chip_id:  # sample_id匹配
                            breed_ped_data.append(sample_data)
                            id_mapping.append(f"{chip_id}\t{sample_name}")
                            found = True
                            break
                    if found:
                        break

            # 写入PED文件
            with open(ped_filename, 'w', encoding='utf-8') as ped_file:
                for sample_data in breed_ped_data:
                    sample_id, family_id, father_id, mother_id, sex, phenotype, genotypes = sample_data
                    line = f"{family_id}\t{sample_id}\t{father_id}\t{mother_id}\t{sex}\t{phenotype}"
                    for genotype in genotypes:
                        line += f"\t{genotype}"
                    ped_file.write(line + "\n")

            # 写入ID映射文件
            with open(idmap_filename, 'w', encoding='utf-8') as idmap_file:
                for mapping in id_mapping:
                    idmap_file.write(mapping + "\n")

            # 创建压缩包
            with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
                zipf.write(ped_filename)
                zipf.write(map_filename)
                zipf.write(idmap_filename)

            # 删除临时文件
            os.remove(ped_filename)
            os.remove(map_filename)
            os.remove(idmap_filename)

            print(f"已生成品种 {breed} 的数据文件: {zip_filename}")

        except Exception as e:
            print(f"生成品种 {breed} 文件时出错: {e}")
            self.error_status = True

    def run(self):
        """运行主程序"""
        print("程序开始运行，请稍候")

        # 1. 检查软件有效期（注释掉以便测试）
        # self.check_expiry()

        # 2. 获取用户输入
        while True:
            farm_code = input("请输入四位场代码: ").strip()
            if len(farm_code) == 4:
                break
            print("输入的场代码不是4位，请重新输入")

        # 3. 加载配置文件
        self.load_exclude_list()
        self.load_split_table()

        if self.error_status:
            print("请检查并修正问题后重新运行程序")
            input("请按任意键退出窗口")
            sys.exit(1)

        # 4. 查找并加载PLINK文件
        plink_files = self.find_plink_files()
        if self.error_status:
            print("请检查并修正问题后重新运行程序")
            input("请按任意键退出窗口")
            sys.exit(1)

        # 5. 加载所有MAP和PED文件
        for map_file, ped_file in plink_files:
            prefix = map_file[:-4]

            # 加载MAP文件
            map_data = self.load_map_file(map_file)
            if map_data:
                self.map_data[prefix] = map_data

            # 加载PED文件
            if map_data:
                ped_data = self.load_ped_file(ped_file, len(map_data))
                if ped_data:
                    self.ped_data[prefix] = ped_data

        if self.error_status:
            print("请检查并修正问题后重新运行程序")
            input("请按任意键退出窗口")
            sys.exit(1)

        # 6. 验证数据一致性
        if not self.validate_map_consistency():
            print("请检查并修正问题后重新运行程序")
            input("请按任意键退出窗口")
            sys.exit(1)

        # 7. 验证芯片格式（可选）
        # self.validate_v1plus_chip()

        # 8. 验证数据完整性
        if not self.validate_data_integrity():
            print("请检查并修正问题后重新运行程序")
            input("请按任意键退出窗口")
            sys.exit(1)

        # 9. 按品种拆分数据
        self.split_by_breed()

        if not self.error_status:
            print("恭喜，成功拆分基因型数据！")
        else:
            print("程序执行过程中出现错误")

        input("请按任意键退出窗口")


def main():
    """主函数"""
    splitter = GenotypeDataSplitter()
    splitter.run()


if __name__ == "__main__":
    main()
