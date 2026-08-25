import re
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("file_name", help="Input file name")
args = parser.parse_args()

# 1. 讀取原始資料
with open(args.file_name, 'r', encoding='utf-8') as f:
    lines = f.readlines()

# 2. 準備正則表達式，用來擷取 Gamma, Beta, Tmp Cut 數值
pattern = re.compile(r'Gamma:\s*([0-9.+-Ee]+),\s*Beta:\s*([0-9.+-Ee]+),\s*Max Cut:\s*([0-9.+-Ee]+)')

# 3. 建立一個 list 來存放解析後的資料
parsed_data = []
for line in lines:
    match = pattern.search(line)
    if match:
        gamma_str, beta_str, tmp_cut_str = match.groups()
        parsed_data.append([gamma_str, beta_str, tmp_cut_str])

# 3.1 將字串轉成浮點數後進行排序 (先 Gamma，再 Beta)
#     如果 Gamma 或 Beta 是整數也沒關係，float() 同樣能處理。
parsed_data.sort(key=lambda x: (float(x[0]), float(x[1])))

# 4. 將資料寫入 CSV 檔案 (Excel 可讀取)
with open(args.file_name+'.csv', 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    # 寫入欄位名稱
    writer.writerow(["Gamma", "Beta", "Tmp Cut"])
    # 寫入每一行資料
    writer.writerows(parsed_data)

print(f"解析並排序完成，結果已存為 {args.file_name+'.csv'}。")
