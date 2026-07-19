from MKmodel import SMK, ModifiedSMK

import numpy as np
import pandas as pd
import scipy.optimize as opt
import matplotlib.pyplot as plt
import os

def hlet_survival(D, alpha, beta):
    '''
    LQ 拟合函数
    '''
    return np.exp(-alpha * D - beta * D**2)

def cal_alpha0_beta0(file, iscorrected=True):
    '''
    计算低LET的alpha0和beta0
    '''
    data = pd.read_csv(file)
    if iscorrected:
        dose = data['nominal_dose_Gy'].values * 0.9
    else:
        dose = data['nominal_dose_Gy'].values
    survival = data['survival_fraction'].values
    ln_survival = np.log(survival)
    def lqfit_func(x, alpha, beta):
        return -alpha * x - beta * x**2

    # 使用非线性最小二乘法拟合LQ模型
    popt, pcov = opt.curve_fit(lqfit_func, dose, ln_survival, p0=(0.1, 0.01))
    alpha0, beta0 = popt
    return alpha0, beta0

def loss_function(survival_pred, survival_data):
    '''
    计算损失函数
    '''
    return np.sum((np.log(survival_pred) - np.log(survival_data)) ** 2)

def cal_loss_of_z0(rootfile, branch, dose, r_d, alpha_hlet, beta_hlet, alpha_0, beta_0, z0_set, outdir):
    '''
    计算不同z0下的损失函数
    '''

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    data_carbon = hlet_survival(dose, alpha_hlet, beta_hlet)
    plt.plot(dose, data_carbon, 'o-', label='Carbon Data', color='black')
    loss_values = []
    print("="*100)
    print(f"Calculating loss for r_d = {r_d} nm")
    print("="*100)
    print(f"| {'z0':<5} | {'Loss':<10} | {'z_d_F':<10} | {'z_n_F':<10} | {'z_d_D':<10} | {'z_n_D':<10} | {'zs_d_F':<10} | {'zs_d_D':<10} |")
    for z0 in z0_set:
        # smk = SMK.from_root(rootfile, branch, alpha0=alpha_0, beta0=beta_0, z0=z0)
        # survival_pred = smk.survival(dose, N=100000)
        # data_carbon = hlet_survival(dose, alpha_hlet, beta_hlet)
        # loss = loss_function(survival_pred, data_carbon)
        # loss_values.append(loss)
        # print(f"| {z0:<5} | {loss:<10.3e} | {smk.z_d_F:<10.4f} | {smk.z_n_F:<10.4f} | {smk.z_d_D:<10.4f} | {smk.z_n_D:<10.4f} | {smk.zs_d_F:<10.4f} | {smk.zs_d_D:<10.4f} |")
        # plt.plot(dose, survival_pred, '-',label=f'z0={z0} Gy')

        msmk = ModifiedSMK.from_root(rootfile, branch, alpha0=alpha_0, beta0=beta_0, z0=z0)
        survival_pred = msmk.survival(dose)
        data_carbon = hlet_survival(dose, alpha_hlet, beta_hlet)
        loss = loss_function(survival_pred, data_carbon)
        loss_values.append(loss)
        print(f"| {z0:<5} | {loss:<10.3e} | {msmk.z_d_F:<10.4f} | {msmk.z_n_F:<10.4f} | {msmk.z_d_D:<10.4f} | {msmk.z_n_D:<10.4f} | {msmk.zs_d_D:<10.4f} |")
        plt.plot(dose, survival_pred, '-',label=f'z0={z0} Gy')


    print("="*100)
    plt.yscale('log')
    plt.title(f"Carbon Survival Curve Fit for r_d = {r_d} nm")
    plt.xlabel("Dose (Gy)")
    plt.ylabel("Survival Fraction")
    plt.ylim(min(data_carbon)/2, 1.2)
    plt.legend()
    plt.savefig(os.path.join(outdir, f"carbon_survival_rd_{r_d}nm.png"))
    plt.close()

    return loss_values


def main():
    ### 输入参数
    alpha_hlet = 1.82
    beta_hlet = 0.169
    r_n = 11.2
    r_d = 0.324
    z0_set = np.arange(60, 150, 10)
    low_let_file = "./data/validation/Xray_survival.csv"
    high_let_dir = "./data/Am241"
    # high_let_file = "data/Am241/am241_phy_decay_rd200nm.root"
    out_dir = "./result/carbon_fit"
    tree = "single_events"

    ### 计算低LET的alpha0和beta0
    alpha_0, beta_0 = cal_alpha0_beta0(low_let_file, iscorrected=True)
    print(f"Low LET alpha0: {alpha_0}, beta0: {beta_0}")

    all_top10_per_rd = []

    # 遍历每一个r_d的高LET文件，计算不同z0下的损失函数，并找到最佳的z0值
    for file in sorted(os.listdir(high_let_dir)):
        if file.endswith(".root"):
            high_let_file = os.path.join(high_let_dir, file)
            r_d = float(file.split("_rd")[1].split("nm")[0])

        ### 计算不同z0下的损失函数，以及绘制生存曲线
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        dose = np.linspace(0, 1, 20)
        loss_values = cal_loss_of_z0(high_let_file, tree, dose, r_d, alpha_hlet, beta_hlet, alpha_0, beta_0, z0_set, out_dir)

        # 1. 打包当前rd下所有(z0, loss)配对
        z_loss_pairs = list(zip(loss_values, z0_set))
        # 2. 按loss从小到大排序
        z_loss_sorted = sorted(z_loss_pairs, key=lambda x: x[0])
        # 3. 取当前rd内部loss最小前10
        top10_current_rd = z_loss_sorted[:10]

        print(f"==== r_d = {r_d} nm 内部loss前10 ====")
        for rank, (loss, z0) in enumerate(top10_current_rd, 1):
            # print(f"第{rank:2d} | z0 = {z0}, Loss = {loss:.6f}")
            # 存入总列表：(loss, r_d, z0)
            all_top10_per_rd.append((loss, r_d, z0))
            
    
    # 打包成 (loss, r_d, z0) 元组列表
    data_sorted = sorted(all_top10_per_rd, key=lambda x: x[0])
    top5 = data_sorted[:15]

    # 循环打印前10
    for idx, (loss, rd, z0) in enumerate(top5, 1):
        print(f"第{idx:2d}优: r_d = {rd} nm, z0 = {z0}, Loss = {loss}")

if __name__ == "__main__":
    main()

