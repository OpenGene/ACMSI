import io
import os
import cv2
import pickle
import argparse
import numpy as np
import pandas as pd
from PIL import Image
from skimage import morphology
import matplotlib.pyplot as plt
# from line_profiler import LineProfiler
from skimage.metrics import mean_squared_error as mse


## 原参数：YM：100；NPM：6；EPM：15；EPR：44
## 优化后：YM：30；NPM：7；EPM：40；EPR：40
def preprocess(df,min_val,max_val):
    df = df.query(f"Size >= {min_val} and Size <= {max_val}")
    df = df.reset_index()
    # rows_to_drop = []  # 存储需要删除的行的索引
    # if max_val != 320:
    #     for idx, row in df.iterrows():
    #         if (idx > 1) and (idx < (len(df) - 2)):
    #             # 行位于列表中间，比对上下三行的差值，要求<50倍
    #             if (row['Height'] > 50*df.iloc[idx - 2]['Height']
    #                     and row['Height'] > 50*df.iloc[idx - 1]['Height']
    #                     and row['Height'] > 50*df.iloc[idx + 1]['Height']
    #                     and row['Height'] > 50*df.iloc[idx + 2]['Height']):
    #                 rows_to_drop.append(idx)
    # df = df.drop(rows_to_drop)  # 删除需要删除的行
    # df.loc[df['Height'] > 3000,'Height'] = 3000
    return df

def find_NT(df, Nname, Tname, marker, min_val, max_val):
    global color_list
    N = df[df["sample"]==Nname]
    T = df[df["sample"]==Tname]
    N_color = N[N['channel'].str.contains(marker)]
    T_color = T[T['channel'].str.contains(marker)]
    ## 丢弃异常值
    N_color = preprocess(N_color,min_val,max_val)
    T_color = preprocess(T_color,min_val,max_val)
    return N_color, T_color

### 生成打分图像
def create_image(df, min_val, max_val, k):
    fig = plt.figure(figsize=(4, 3))
    # 绘制柱状图
    plt.bar(df['Size']+k, df['Height'],align='center')
    plt.xlim(min_val, max_val)
    plt.ylim(30)
    plt.axis('off')
    buffer = io.BytesIO()
    plt.savefig(buffer, format='PNG')
    plt.close(fig)  # 关闭图形对象
    buffer.seek(0)
    
    filtered_df = df[(df['Size'] >= min_val) & (df['Size'] <= max_val)] 
    max_area = filtered_df['Height'].max()  # 获取筛选后的数据中最大的area值
    return buffer, max_area

### 作图寻优
def create_images(df, min_val, max_val, k):
    fig = plt.figure(figsize=(4, 3))
    # 绘制柱状图
    plt.bar(df['Size']+k, df['Height'], align='center')
    plt.xlim(min_val, max_val)
    plt.ylim(100)
    plt.axis('off')
    buffer = io.BytesIO()
    plt.savefig(buffer, format='PNG')
    plt.close(fig)  # 关闭图形对象
    buffer.seek(0)

    return buffer

def plot_images_gene(image1, image2, Tname, Nname, gene):
    global color_list
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
    plt.suptitle(f"Gene: {gene}", fontsize=12, x=0.10, y=0.85)
    axes[1].imshow(image1, cmap='gray')
    axes[0].imshow(image2, cmap='gray')
    # 取消坐标轴
    axes[0].axis('off')
    axes[1].axis('off')
    axes[0].set_title(f"Tumor: {Tname}")
    axes[1].set_title(f"Normal: {Nname}")
    plt.show()

def show_array(array):
    # with open("real-MSIresult.csv", "a") as file:
    with open("image_array.csv", "a") as file:
        for row in array:
            row_str = ' '.join(str(x) for x in row)
            file.write("KKKKKK"+ row_str + '\n')

def plot_sub_img(image1, image2, Tname, gene, ax, best_k, num_peaks, num_height, max_height):
    score, num_height, num_peaks, peak_image, height_image = msi_score(image1,image2)
    new_image = 2* peak_image + image2 + height_image
    ax.imshow(new_image, cmap = 'magma', vmin=552, vmax=1020)  # 使用 'viridis' colormap，您可以根据需要选择不同的colormap
    ax.set_title(f"{gene}")
    ax.text(10, 70, f"Max peaks: {max_height:.0f} RFU", color='black', fontsize=10)
    ax.text(10, 30, f"New peaks: {num_peaks}", color='black', fontsize=10)
    ax.text(10, 50, f"Extra peaks: {num_height}", color='black', fontsize=10)
    ax.text(10, 10, f"Shift: {best_k/7:.2f} bp", color='black', fontsize=10)
    # ax.colorbar()  # 添加颜色条
    ax.axis('off')

def msi_score(image1,image2):
    # 定义结构元素
    
    selem_h = morphology.rectangle(6, 3)
    selem_hp = morphology.rectangle(40, 3)
    selem_np = morphology.rectangle(7, 3)


    # 获取每组图片的非白色像素数量[0,100,156(image1-image2)]
    dif_hp = image2 - image1
    dif_hp = np.where(dif_hp == 100, 99, 0)
    dif_np = image2 - image1
    dif_np = np.where(dif_np == 100, 99, 0)

    ##  extra peaks
    bott = dif_hp[-40:-34, :]
    zero_cols_hp = np.all(bott == 0, axis=0)
    dif_hp[:, zero_cols_hp] = 0
    # 对于extra peaks要求严格，必须增长到一定程度，呈增长趋势
    peaks_hp = morphology.opening(dif_hp, selem_hp)

    ## new peaks
    zero_cols_np = np.where(dif_np[-34, :] == 0)[0]
    # 将这些列的所有元素设为0
    dif_np[:, zero_cols_np] = 0
    peaks_np = morphology.opening(dif_np, selem_np)

    peaks_dif =  peaks_np + peaks_hp
    peak_image = np.where(peaks_dif == 0, 255, 99)

    # height_dif_tmp = np.clip((image2 - image1), 0, None)
    height_dif_tmp = image2 -image1
    height_dif_tmp = np.where(height_dif_tmp == 100, 99, 255)
    height_dif =  abs(height_dif_tmp - peak_image)

    height_dif = morphology.opening(height_dif, selem_h)
    height_image = np.where(height_dif == 156, 99, 255)
    

    num_peaks = np.sum(peak_image != 255)
    num_height = np.sum(height_image != 255)
    score = 0.7*num_peaks + 0.3*num_height
    return score, num_height, num_peaks, peak_image, height_image

## 用整体align
def fitting_median(sampleN, sampleT, min_val, max_val, rp):
    buffer2 = create_images(sampleT, min_val, max_val, 0)
    image_data2 = np.asarray(bytearray(buffer2.read()), dtype=np.uint8)
    image2 = cv2.imdecode(image_data2, cv2.IMREAD_GRAYSCALE)
    
    left = -rp
    right = rp
    
    while right - left > 0.01:
        mid1 = left + (right - left) / 3
        mid2 = right - (right - left) / 3
        
        buffer_mid1 = create_images(sampleN, min_val, max_val, mid1)
        buffer_mid2 = create_images(sampleN, min_val, max_val, mid2)
        
        image_data_mid1 = np.asarray(bytearray(buffer_mid1.read()), dtype=np.uint8)
        image_data_mid2 = np.asarray(bytearray(buffer_mid2.read()), dtype=np.uint8)
        image_mid1 = cv2.imdecode(image_data_mid1, cv2.IMREAD_GRAYSCALE)
        image_mid2 = cv2.imdecode(image_data_mid2, cv2.IMREAD_GRAYSCALE)
        
        mse_mid1 = mse(image_mid1, image2)
        mse_mid2 = mse(image_mid2, image2)
        
        if mse_mid1 < mse_mid2:
            right = mid2
        else:
            left = mid1

    best_k = (left + right) / 2
    best_mse = mse(image_mid1, image2)

    buffer2.close()
    buffer_mid1.close()
    buffer_mid2.close()
    
    return best_k, best_mse

## 对图片直接进行align
def align_image(image1, image2, rp):
    # 定义一个差异度阈值，以便在何时认为两幅图像对齐
    threshold = 0.001  # 可以根据实际情况调整
    ## 根据图片色素计算而来
    step = 7*(2*rp-1)
    # 初始化最小差异度和最佳偏移值
    min_difference = float('inf')
    best_offset = 0
    # 遍历不同的偏移值，进行左右移动
    for offset in range(-step, step):  # 范围可以根据图像大小和移动范围调整
        # 创建一个偏移后的图像
        shifted_image = np.roll(image1, offset, axis=1)

        # 计算差异度
        difference = np.sum(np.abs(image2 - shifted_image))

        # 如果差异度低于阈值，更新最小差异度和最佳偏移值
        if difference < min_difference:
            min_difference = difference
            best_offset = offset

    # 根据最佳偏移值进行图像对齐
    aligned_image = np.roll(image1, best_offset, axis=1)
    return aligned_image, image2, best_offset

def search_gene(color, rp):
    search_table = {
        ('NED', 1): 'BAT25',
        ('NED', 2): 'D17S250',
        ('VIC', 1): 'BAT26',
        ('VIC', 2): 'D2S123',
        ('6-FAM', 2): 'D5S346',
        ('6-FAM', 4): 'pentaC'
    }
    result = search_table.get((color, rp))
    return result

def classify_image(df, Sample, Nname, Tname, site, marker, min_val, max_val, rp, ax):
    global fu_list
    sampleN, sampleT = find_NT(df, Nname, Tname, marker, min_val, max_val)
    if (len(sampleN) <= 20) or (len(sampleT) <= 20) and (abs(len(sampleN) - len(sampleT)) > 20):
        return
    ## 生成打分图像
    buffer1, max_height1 = create_image(sampleN, min_val, max_val, 0)
    image_data1 = np.asarray(bytearray(buffer1.read()), dtype=np.uint8)
    image1_orign = cv2.imdecode(image_data1, cv2.IMREAD_GRAYSCALE)
    buffer2, max_height2 = create_image(sampleT, min_val, max_val, 0)
    image_data2 = np.asarray(bytearray(buffer2.read()), dtype=np.uint8)
    image2_orign = cv2.imdecode(image_data2, cv2.IMREAD_GRAYSCALE)
    ## 比对
    image1, image2, best_k = align_image(image1_orign, image2_orign, rp)

    # 下部分切割
    image1_lower = image1[image1.shape[0]//2:, :]
    image2_lower = image2[image2.shape[0]//2:, :]

    best_mse = mse(image1, image2)
    score, num_height, num_peaks, peak_image, height_image = msi_score(image1_lower, image2_lower)
    max_height = min(max_height1, max_height2)

    fu_list.append([f"{Sample} ({Tname}-{Nname})",site, round(max_height, 0), round(best_k/7, 2),num_peaks])
    print(f"Sample: {Nname}\t Gene: {site}\nK: {best_k}, MSE: {best_mse}, New Peaks: {num_peaks}, New Height: {num_height}, Max Height: {max_height}, Score: {score}\n")
    plot_sub_img(image1, image2, Tname, site, ax, best_k, num_peaks, num_height, max_height)

def main():
    # 设置工作目录为当前脚本所在目录
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    global fu_list
    with open('www/MSI_SVM.pkl', 'rb') as f:
        svm_model = pickle.load(f)

    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--TP", type=argparse.FileType('r'), nargs='?', help="input TP-information file")
    parser.add_argument("-p","--peak", type=argparse.FileType('r'), nargs='?', help="input Peaks-information file")
    parser.add_argument("-a","--anno", type=argparse.FileType('r'), nargs='?', help="input Annotation-information file")
    args = parser.parse_args()

    TPinfo = args.TP
    peak = args.peak
    anno = args.anno
    df = pd.read_csv(peak, sep='\t', names=['sample','channel','Size','Height'])
    anno = pd.read_csv(anno)
    info = pd.read_csv(TPinfo)
    df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    sampleNum = len(info)
    fig_all, ax_all = plt.subplots(sampleNum , 1,figsize=(18,3*sampleNum ))

    for idx1, row1 in info.iterrows():
        fig, axes = plt.subplots(1, 6, figsize=(18, 3))
        Nname = row1['Nname']
        Tname = row1['Tname']
        Sample = row1['Sample']
        for idx2, row2 in anno.iterrows():
            site = row2['site']
            marker = row2['marker']
            range_str = row2['range']
            range_list = [float(x) for x in range_str.split('-')]
            size = row2['size']
            classify_image(df, Sample, Nname, Tname, site, marker, range_list[0], range_list[1], size, axes[idx2])
        fig.suptitle(f"{Sample} ({Tname}-{Nname})")
        plt.tight_layout(pad=1.5)

        buf = io.BytesIO()
        fig.savefig(buf, format='png')
        buf.seek(0)

        # 从缓存中读取子图加载到大图中
        img = Image.open(buf)
        ax_all[idx1].imshow(img)
        ax_all[idx1].axis('off')

        buf.close() 
        plt.close(fig) 
    plt.savefig("www/plot_all_img.png")
    plt.close()
    fu_df = pd.DataFrame(fu_list, columns=['sample','site', 'RFU', 'shift', 'peaks'])
    # 进行预测
    X = fu_df['peaks'].values.reshape(-1, 1)
    predictions = svm_model.predict(X)
    fu_df['mark'] = predictions
    fu_df.to_csv('www/quality_info.csv', index=False)

if __name__ =='__main__':
    fu_list = []
    main()
    # lp = LineProfiler()
    # lp.add_function(preprocess)
    # lp.add_function(find_NT)
    # lp.add_function(create_image)
    # lp.add_function(create_images)
    # lp.add_function(plot_images_gene)
    # lp.add_function(plot_sub_img)
    # lp.add_function(msi_score)
    # lp.add_function(fitting_median)
    # lp.add_function(align_image)
    # lp.add_function(search_gene)
    # lp.add_function(classify_image)

    # test_func = lp(main)
    # test_func()
    # lp.print_stats()