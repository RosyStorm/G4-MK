from MKmodel import SMK

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

def cal_loss_of_z0(rootfile, branch, dose, alpha_hlet, beta_hlet, alpha_0, beta_0, z0_set, outdir):
    '''
    计算不同z0下的损失函数
    '''
    data_carbon = hlet_survival(dose, alpha_hlet, beta_hlet)
    plt.plot(dose, data_carbon, 'o', label='Carbon Data', color='black')
    loss_values = []
    print(f"| {'z0':<5} | {'Loss':<10} | {'z_d_F':<10} | {'z_n_F':<10} | {'z_d_D':<10} | {'z_n_D':<10} | {'zs_d_F':<10} | {'zs_d_D':<10} |")
    for z0 in z0_set:
        smk = SMK.from_root(rootfile, branch, alpha0=alpha_0, beta0=beta_0, z0=z0)
        survival_pred = smk.survival(dose, N=10000)
        data_carbon = hlet_survival(dose, alpha_hlet, beta_hlet)
        loss = loss_function(survival_pred, data_carbon)
        loss_values.append(loss)
        print(f"| {z0:<5} | {loss:<10.3e} | {smk.z_d_F:<10.4f} | {smk.z_n_F:<10.4f} | {smk.z_d_D:<10.4f} | {smk.z_n_D:<10.4f} | {smk.zs_d_F:<10.4f} | {smk.zs_d_D:<10.4f} |")
        plt.plot(dose, survival_pred, '-',label=f'z0={z0} Gy')

    plt.yscale('log')
    # plt.legend()
    plt.savefig(os.path.join(outdir, "carbon_survival.png"))

    return loss_values


def main():
    ### 输入参数
    alpha_hlet = 1.82
    beta_hlet = 0.169
    r_n = 11.2
    r_d = 0.324
    z0_set = np.arange(100, 120, 2)
    low_let_file = "./data/validation/Xray_survival.csv"
    high_let_file = "data/am241_phy_decay_Membrane.root"
    out_dir = "./result/carbon_fit"
    tree = "single_events"

    ### 计算低LET的alpha0和beta0
    alpha_0, beta_0 = cal_alpha0_beta0(low_let_file, iscorrected=True)
    print(f"Low LET alpha0: {alpha_0}, beta0: {beta_0}")

    ### 计算不同z0下的损失函数，以及绘制生存曲线
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    dose = np.linspace(0, 1, 20)
    loss_values = cal_loss_of_z0(high_let_file, tree, dose, alpha_hlet, beta_hlet, alpha_0, beta_0, z0_set, out_dir)

    ### 找到最佳的z0值
    best_loss = min(loss_values)
    best_z0 = z0_set[loss_values.index(best_loss)]
    print(f"Best z0: {best_z0}, Best Loss: {best_loss}")


if __name__ == "__main__":
    main()

