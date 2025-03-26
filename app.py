from flask import Flask, send_file, request, render_template, redirect, url_for, session
import os
import io
from werkzeug.utils import secure_filename
import pandas as pd
from datetime import datetime
import matplotlib
import matplotlib.pyplot as plt
import base64

app = Flask(__name__)
matplotlib.use('Agg')  # 使用非交互式后端

# 配置上传文件夹和允许的文件类型
UPLOAD_FOLDER = './'
ALLOWED_EXTENSIONS = {'txt', 'csv'}

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
# 必须设置 secret_key 才能使用 session
app.secret_key = 'mpquic'
def calculate_fold_change(target_gene, processed_gene_list, internal_reference, group, control_sample):
    for gene_name in target_gene:
        name = gene_name + '/' + internal_reference
        if name not in processed_gene_list:
            processed_gene_list.append(name)
        control_value = group.loc[group['Sample Name'] == control_sample, name].values[0]
        # 计算Fold Change
        changed_name = 'Flod Change ' + name
        if changed_name not in processed_gene_list:
            processed_gene_list.append(changed_name)
        group[changed_name] = group[name] / control_value
    return group
def allowed_file(filename):
    """检查文件扩展名是否允许"""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')

@app.route('/upload', methods=['GET', 'POST'])
def upload_files():
    if request.method == 'POST':
        if 'files' not in request.files:
            return "Error: No file part in the request", 400

        files = request.files.getlist('files')  # 获取多个文件
        saved_files = []

        for file in files:
            if file.filename == '':
                continue  # 忽略空文件

            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                saved_files.append(filename)
            else:
                return f"Error: Invalid file type ({file.filename}). Only TXT and CSV are allowed.", 400

        if saved_files:
            print(f'file list: {saved_files}')
            # 使用列表推导式读取所有文件
            dfs = [pd.read_csv(file_name, sep="\t") for file_name in saved_files]
            # 合并所有数据
            df = pd.concat(dfs, ignore_index=True)
            # 删除不需要的行
            df = df[df['Sample Name'].str.match(r'^[0-9]')]
            # 删除 Concentration Mean 列中值为 '-' 的行
            df = df[df['Concentration Mean'] != "-"]
            # 保留Sample Name和Concentration Mean数据
            df = df.pivot_table(index='Sample Name', columns='Gene Name', values='Concentration Mean')
            # 重置索引，确保Sample Name成为列而不是索引
            df.reset_index(inplace=True)

            # 将 DataFrame 转换为 CSV 格式存储在 session 中
            csv_data = df.to_csv(index=False)
            session['df'] = csv_data  # 存储在 session 中

            for del_file in saved_files:
                # 文件处理后立即删除
                os.remove(del_file)  # 提取数据后删除文件

            return redirect(url_for('input_gene'))
        else:
            return "Error: No valid files uploaded", 400

    return render_template('upload.html')

@app.route('/input_gene', methods=['GET', 'POST'])
def input_gene():
    # 从 session 获取 CSV 数据
    csv_data = session.get('df')
    if csv_data:
        # 将 CSV 数据转换回 DataFrame
        df = pd.read_csv(io.StringIO(csv_data))
        # 获取表头 (列名)
        column_names = df.columns.tolist()
        if 'Sample Name' in column_names:
            column_names.remove('Sample Name')

        if request.method == 'POST':
            # 获取用户输入
            internal_reference = request.form['internal_reference']
            # 只将非空的 Target Gene 输入框的值添加到列表
            target_genes = [
                request.form.get(f'target_gene_{i}')
                for i in range(1, 6)
                if request.form.get(f'target_gene_{i}')
            ]
            ori_gene_list = target_genes + [internal_reference]

            # 打印用户输入的数据，或做进一步的处理
            # print("Original gene list:", ori_gene_list)
            print("Internal Reference:", internal_reference)
            print("Target Genes:", target_genes)

            for i in target_genes:
                # 创建新列 target_gene/internal_gene
                df[i + '/' + internal_reference] = df[i] / df[internal_reference]
            # 添加Round列，根据Sample Name中的数字分配Round
            df['Round'] = df['Sample Name'].str.extract(r'(\d)')[0].astype(int)
            # 修改Sample Name列：去掉数字前缀并替换指定值
            df['Sample Name'] = df['Sample Name'].str.replace(r'(^\d)', '', regex=True)  # 只去除数字前缀

            # 将 DataFrame 转换为 CSV 格式存储在 session 中
            csv_data = df.to_csv(index=False)
            session['df'] = csv_data  # 存储在 session 中
            # 将target_genes ori_gene列表转换为逗号分隔的字符串
            target_genes_str = ','.join(target_genes)
            ori_gene_list_str = ','.join(ori_gene_list)
            return redirect(url_for('rename_samples', target_genes=target_genes_str, ori_gene_list=ori_gene_list_str, internal_reference=internal_reference))
    else:
        return "No data available", 400  # 如果 session 中没有数据，返回错误

    # 显示表头和输入框的页面
    return render_template('input_gene.html', column_names=column_names)

@app.route('/rename_samples', methods=['GET', 'POST'])
def rename_samples():
    # 从input_gene获取变量
    internal_reference = request.args.get('internal_reference')
    ori_gene_list_str = request.args.get('ori_gene_list')
    ori_gene_list = ori_gene_list_str.split(',')
    target_genes_str = request.args.get('target_genes')
    target_genes = target_genes_str.split(',')
    # 从 session 获取 CSV 数据
    csv_data = session.get('df')
    if csv_data:
        # 将 CSV 数据转换回 DataFrame
        df = pd.read_csv(io.StringIO(csv_data))

        # 获取唯一的旧名称列表
        sample_names = df['Sample Name'].unique()

        if request.method == 'POST':
            rename_dict = {}

            # 获取用户输入的新名称并生成字典
            for name in sample_names:
                new_name = request.form.get(f'rename_{name}')
                if new_name:  # 确保用户输入了新的名字
                    rename_dict[name] = new_name

            # 获取用户选择的 base_name
            base_name = request.form.get('base_name')
            base_name_value = rename_dict.get(base_name)
            print(f"Selected base name: {base_name_value}")
            if base_name:
                # 以 base_name 为参考，更新其他 Sample Name
                base_name_value = rename_dict.get(base_name)
                if base_name_value:
                    for name, new_name in rename_dict.items():
                        if name != base_name:
                            # 对其他项进行处理，使用 base_name 作为参考（例如，给它们添加后缀等）
                            rename_dict[name] = new_name
            # print(f"rename list: {rename_dict}")
            # 批量替换
            df['Sample Name'] = df['Sample Name'].replace(rename_dict)
            # 对每个Round进行分组并计算Fold Change
            processed_gene_list = []
            # print(f"target gene:", target_genes)
            # print(f"internal_reference", internal_reference)
            # print(f"base_name_value", base_name_value)
            # print(f"oringnal gene list:", ori_gene_list)
            df = df.groupby('Round').apply(
                lambda x: calculate_fold_change(target_genes, processed_gene_list, internal_reference, x, base_name_value))

            # 确保 Round 列是普通列，而不是索引的一部分
            df['Round'] = df.index.get_level_values('Round')

            # 重置索引，将 Round 列移回到数据列
            df.reset_index(drop=True, inplace=True)

            # 只保留四列，确保Fold Change列存在
            for i in ori_gene_list:
                df = df.drop(columns=[i])

            # 将数据分为四个表
            df_dict = {}
            for i in processed_gene_list:
                df_dict[i] = df.pivot_table(index='Sample Name', columns='Round', values=i)
                df_dict[i]['Average'] = df_dict[i].iloc[:, 1:].mean(axis=1)

            # 生成柱状图并转换为 Base64
            plt.style.use('seaborn-darkgrid')  # 选择冷淡风格
            img_data_list = []
            for gene, df_table in df_dict.items():
                fig, ax = plt.subplots(figsize=(8, 5))

                # 转置 DataFrame，使 Round 成为 x 轴，Sample Name 作为不同类别
                df_table.T.iloc[:-1].plot(kind='bar', ax=ax, colormap='Blues', edgecolor='black')  # 去掉 'Average' 列
                plt.title(f"{gene}")
                plt.xlabel("Round")
                plt.ylabel("Expression Level")
                plt.xticks(rotation=0)  # x 轴标签水平显示
                plt.legend()

                # 将图例放置在表格外右侧上方
                ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), frameon=False)

                plt.tight_layout()  # 自动调整布局，避免图例超出边界

                # 保存图片到内存
                img_io = io.BytesIO()
                plt.savefig(img_io, format='png', bbox_inches='tight')
                img_io.seek(0)
                img_data = base64.b64encode(img_io.getvalue()).decode('utf-8')
                img_data_list.append(img_data)
                plt.close(fig)  # 关闭图像，防止内存占用

            # 删除当前文件夹下所有csv文件 获取当前工作目录
            current_directory = os.getcwd()
            # 遍历当前目录下的所有文件
            for filename in os.listdir(current_directory):
                if filename.endswith('.csv'):
                    file_path = os.path.join(current_directory, filename)
                    try:
                        # 删除文件
                        os.remove(file_path)
                        print(f"Deleted: {file_path}")
                    except Exception as e:
                        print(f"Error deleting file {file_path}: {e}")

            # 获取当前日期
            current_date = datetime.now().strftime("%Y%m%d")
            # 生成文件名
            file_name = f"{current_date}_{'_'.join(ori_gene_list)}.csv"
            # file_name = f"output.csv"

            # 创建并写入 CSV 文件
            with open(file_name, mode='w', newline='') as f:
                for i in df_dict:
                    # 添加空行分隔
                    f.write("\n")
                    f.write(f"{i}\n")
                    df_dict[i].to_csv(f, header=True)

            # 读取 CSV 文件
            df_csv = pd.read_csv(file_name)
            # 将 DataFrame 转换为 HTML 表格
            table_html = df_csv.to_html(classes='table table-striped', index=False)

            # 返回 CSV 内容并在前端显示，同时提供下载
            # return render_template('download.html', table_html=table_html, file_name=file_name)
            return render_template('download.html', file_name=file_name, img_data_list=img_data_list)
    else:
        return "No data available", 400  # 如果 session 中没有数据，返回错误

    return render_template('rename_samples.html', sample_names=sample_names)

@app.route('/download/<file_name>')
def download(file_name):
    file_path = os.path.join(UPLOAD_FOLDER, file_name)
    # 检查文件是否存在
    if not os.path.exists(file_path):
        return render_template('error.html', error_message="ERROR: File not found!")
    # 返回文件供用户下载
    return send_file(file_name, as_attachment=True, download_name=file_name)

if __name__ == '__main__':
    # 确保上传目录存在
    os.makedirs(UPLOAD_FOLDER, exist_ok=True)
    # app.run(
    #     host='0.0.0.0',
    #     port=8888,  # 默认的 HTTPS 端口
    #     debug=True,
    #     ssl_context=('./ssl/server.crt', './ssl/server.key')  # 配置证书和私钥
    # )
    app.run(
        host='0.0.0.0',
        port=5000,  # Flask 默认 HTTP 端口
        debug=True
    )
