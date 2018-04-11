#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import ode, odeint
import matplotlib.pyplot as plt

u"K-T model(野本の1次系近似モデル:T2=T3=0)のシミュレータ"


def maneuver_u_fix(K, T, x0, y0, psi0, u0, r0, duration, sampling, delta_list):
    u"""
    KT modelでの操縦シミュレーション
    :param K: 操縦性パラメータ[1/s]
    :param T: 操縦性パラメータ[m]
    :param x0:船体の初期位置(x座標[m])
    :param y0:船体の初期位置(y座標[m])
    :param psi0: 船体の初期の向き[rad]
    :param u0: 船体座標系における船体の初期速度[m/s]
    :param r0: 船体座標系における船体の角速度[rad/s]
    :param duration : 操縦シミュレーションの期間[s]
    :param sampling : 時間の区切り個数[-]
    :param delta_list : 舵角情報のリスト(総時間はduration, 登録されている個数はsamplingのlist)
    :return time_list
    :return X
    """
    time_list = np.linspace(0.00, duration, sampling)
    X_init = np.array([x0, y0, psi0, u0, r0])

    def __eom(X, t):
        u"""
        Equation of motion
        :param X:
        :param t:
        :return: [d_x, d_y, d_psi, d_u, d_r]
        """
        # tに対応した尤も確からしい舵角を算出する
        # tに最も近いtime_list内の情報を探し、その時間における舵角情報を採用する
        nearest_id = np.abs(np.asarray(time_list) - t).argmin()
        delta_t = delta_list[nearest_id]

        # 常微分方程式
        d_x = X[3] * np.cos(X[2])  # dx = U*cos(psi)
        d_y = X[3] * np.sin(X[2])  # dy = U*sin(psi)
        d_psi = X[4]  # dpsi = r
        d_u = 0.0  # du = 0 (Fixed)
        d_r = 1.0 / T * (- X[4] + K * (delta_t))  # dr = 1/T * ( -r + K * delta_t)
        return [d_x, d_y, d_psi, d_u, d_r]

    X = odeint(__eom, X_init, time_list)
    X = np.c_[delta_list, X]  # deltaの値もXに入れておく(前に結合)

    return time_list, X


def maneuver(K, T, x0, y0, psi0, u_list, r0, duration, sampling, delta_list):
    u"""
    KT modelでの操縦シミュレーションで、uが変化するかつ既知としたもの
    :param K: 操縦性パラメータ[1/s]
    :param T: 操縦性パラメータ[m]
    :param x0:船体の初期位置(x座標[m])
    :param y0:船体の初期位置(y座標[m])
    :param psi0: 船体の初期の向き[rad]
    :param u_list: 船体の速度の時系列情報(list)[m/s]
    :param r0: 船体座標系における船体の角速度[rad/s]
    :param duration : 操縦シミュレーションの期間[s]
    :param sampling : 時間の区切り個数[-]
    :param delta_list : 舵角情報のリスト(総時間はduration, 登録されている個数はsamplingのlist)
    :return time_list
    :return X
    """

    time_list = np.linspace(0.00, duration, sampling)
    X_init = np.array([x0, y0, psi0, r0])

    def __eom(X, t):
        u"""
        Equation of motion
        :param X:
        :param t:
        :return: [d_x, d_y, d_psi, d_r]
        """
        # tに対応した尤も確からしい速度と舵角を算出する
        # tに最も近いtime_list内の情報を探し、その時間における速度情報と舵角情報を採用する
        nearest_id = np.abs(np.asarray(time_list) - t).argmin()
        delta_t = delta_list[nearest_id]
        u_t = u_list[nearest_id]

        # 常微分方程式
        d_x = u_t * np.cos(X[2])  # dx = U*cos(psi)
        d_y = u_t * np.sin(X[2])  # dy = U*sin(psi)
        d_psi = X[3]  # dpsi = r
        d_r = 1.0 / T * (- X[3] + K * (delta_t))  # dr = 1/T * ( -r + K * delta_t)
        return [d_x, d_y, d_psi, d_r]

    _X = odeint(__eom, X_init, time_list)

    # 順番を変更する
    xX = _X.T[0]
    yX = _X.T[1]
    psiX = _X.T[2]
    rX = _X.T[3]

    X = np.c_[delta_list, xX, yX, psiX, u_list, rX]  # deltaの値もXに入れておく(前に結合)

    return time_list, X


def maneuver_self_organization(K_list, T_list, x0, y0, psi0, u_list, r0, duration, sampling, delta_list):
    u"""
    KT modelでの操縦シミュレーションで、uが変化するかつ既知とし、かつKとTが時刻ごとに変化することを想定したもの
    :param K_list: 操縦性パラメータ[1/s]
    :param T_list: 操縦性パラメータ[m]
    :param x0:船体の初期位置(x座標[m])
    :param y0:船体の初期位置(y座標[m])
    :param psi0: 船体の初期の向き[rad]
    :param u_list: 船体の速度の時系列情報(list)[m/s]
    :param r0: 船体座標系における船体の角速度[rad/s]
    :param duration : 操縦シミュレーションの期間[s]
    :param sampling : 時間の区切り個数[-]
    :param delta_list : 舵角情報のリスト(総時間はduration, 登録されている個数はsamplingのlist)
    :return time_list
    :return X
    """

    time_list = np.linspace(0.00, duration, sampling)
    X_init = np.array([x0, y0, psi0, r0])

    def __eom(X, t):
        u"""
        Equation of motion
        :param X:
        :param t:
        :return: [d_x, d_y, d_psi, d_r]
        """
        # tに対応した尤も確からしい速度と舵角を算出する
        # tに最も近いtime_list内の情報を探し、その時間における速度情報と舵角情報を採用する
        nearest_id = np.abs(np.asarray(time_list) - t).argmin()
        delta_t = delta_list[nearest_id]
        u_t = u_list[nearest_id]
        K_t = K_list[nearest_id]
        T_t = T_list[nearest_id]

        # 常微分方程式
        d_x = u_t * np.cos(X[2])  # dx = U*cos(psi)
        d_y = u_t * np.sin(X[2])  # dy = U*sin(psi)
        d_psi = X[3]  # dpsi = r
        d_r = 1.0 / T_t * (- X[3] + K_t * (delta_t))  # dr = 1/T * ( -r + K * delta_t)
        return [d_x, d_y, d_psi, d_r]

    _X = odeint(__eom, X_init, time_list)

    # 順番を変更する
    xX = _X.T[0]
    yX = _X.T[1]
    psiX = _X.T[2]
    rX = _X.T[3]

    X = np.c_[delta_list, xX, yX, psiX, u_list, rX]  # deltaの値もXに入れておく(前に結合)

    return time_list, X


def maneuver_zigzag(K, T, x0, y0, psi0, u0, r0, duration, sampling, delta_a, psi_a, t1):
    u"""
    KT modelでのzigzag操縦シミュレーション(delta_a/psi_a zigzag maneuvering)
    :param K: 操縦性パラメータ[1/s]
    :param T: 操縦性パラメータ[m]
    :param x0:船体の初期位置(x座標[m])
    :param y0:船体の初期位置(y座標[m])
    :param psi0: 船体の初期の向き[rad]
    :param u0: 船体座標系における船体の初期速度[m/s]
    :param r0: 船体座標系における船体の角速度[rad/s]
    :param duration : 操縦シミュレーションの期間[s]
    :param sampling : 時間の区切り個数[-]
    :param delta_a: zigzag運動における舵角[rad]
    :param psi_a: zigzag運動における転舵するべき方位の絶対値[rad]
    :param t1: 舵角を0[rad]からdelta_a[rad]まで変更するときに必要な時間[s](すなわち転舵速度を表す)
    :return: time_list
    :return: X
    """
    time_list = np.linspace(0.00, duration, sampling)
    delta0 = 0.0
    X_init = np.array([delta0, x0, y0, psi0, u0, r0])

    def __eom(X, t):
        u"""
        Equation of motion
        :param X:
        :param t:
        :return: [d_delta, d_x, d_y, d_psi, d_u, d_r]
        """
        # 現状とtに対応したzigzag試験に対応する舵角を算出する
        d_delta = 0.0
        if delta_a >= 0:
            if t <= t1:  # 0 - t1
                d_delta = delta_a / t1
            elif np.fabs(X[3]) <= psi_a:  # t1 - t2
                d_delta = 0.0
            else:
                if X[3] > 0.0:
                    if X[0] > - delta_a:  # t2 - t3 & t6 - t7
                        d_delta = - delta_a / t1
                    else:  # t3 - t4
                        d_delta = 0.0
                else:
                    if X[0] < delta_a:  # t4 - t5
                        d_delta = delta_a / t1
                    else:  # t5 - t6
                        d_delta = 0.0
        else:
            if t <= t1:  # 0 - t1
                d_delta = delta_a / t1
            elif np.fabs(X[3]) <= psi_a:  # t1 - t2
                d_delta = 0.0
            else:
                if X[3] < 0.0:  # t2 - t3 & t6 - t7
                    if X[0] < - delta_a:
                        d_delta = - delta_a / t1
                    else:  # t3 - t4
                        d_delta = 0.0
                else:
                    if X[0] > delta_a:  # t4 - t5
                        d_delta = delta_a / t1
                    else:  # t5 - t6
                        d_delta = 0.0

        # 常微分方程式
        d_x = X[4] * np.cos(X[3])  # dx = U*cos(psi)
        d_y = X[4] * np.sin(X[3])  # dy = U*sin(psi)
        d_psi = X[5]  # dpsi = r
        d_u = 0.0  # du = 0 (Fixed)
        d_r = 1.0 / T * (- X[5] + K * X[0])  # dr = 1/T * ( -r + K * delta_t)
        return [d_delta, d_x, d_y, d_psi, d_u, d_r]

    X = odeint(__eom, X_init, time_list)

    return time_list, X


def draw_trajectory(p, save_file_name):
    u"""
    軌跡図を出力する
    :param p: maneuverメソッドで出力される操縦性シミュレーションの結果のxy座標情報
    :param save_file_name: 軌跡図の保存path
    :return:
    """
    plt.plot(p[0], p[1])
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.savefig(save_file_name)
    plt.close()


def draw_parameter_graph(time_list, delta_list, psi_list, u_list, r_list, save_file_name):
    u"""
    シミュレーション結果の時系列データを図で出力する
    :param time_list: maneuverメソッドで出力される操縦性シミュレーションの結果に対応する時間リスト
    :param delta_list: maneuverメソッドで出力される操縦性シミュレーションの結果の舵角情報リスト
    :param psi_list: maneuverメソッドで出力される操縦性シミュレーションの結果の方位情報リスト
    :param u_list: maneuverメソッドで出力される操縦性シミュレーションの結果の速度情報リスト
    :param r_list: maneuverメソッドで出力される操縦性シミュレーションの結果の角速度情報リスト
    :param save_file_name: 軌跡図の保存path
    :return:
    """
    plt.plot(time_list, delta_list * 180 / np.pi, label='delta [degree]')
    plt.plot(time_list, psi_list * 180 / np.pi, label='psi [degree]')
    plt.plot(time_list, u_list, label='u [m/s]')
    plt.plot(time_list, r_list * 180 / np.pi, label='r [degree/s]')
    plt.legend()
    plt.xlabel('Time[s]')
    plt.savefig(save_file_name)
    plt.close()


if __name__ == '__main__':
    # for test

    #################################
    # parameter
    #################################
    K = 0.155  # [1/s]
    T = 80.5  # [s]
    u0 = 20 * (1852.0 / 3600)  # [m/s] (knot * 1852/3600)
    x0 = 0.0  # [m]
    y0 = 0.0  # [m]
    psi0 = 0.0 * np.pi / 180  # [rad] (degree * pi / 180)
    r0 = 0.0 * np.pi / 180  # [rad/s] ([degree/s] * pi / 180)
    duration = 500  # [s]
    sampling = 1000

    # Test for fixed rudder
    # delta_list = np.ones(sampling) * 5.0 * np.pi / 180 # [rad]

    # Test for trigonometric rudder
    time_list = np.linspace(0.00, duration, sampling)
    delta_list = np.sin(2.0 * np.pi / 50 * time_list) * 5.0 * np.pi / 180  # [rad](delta0=10.0[degree], Ts=50[s])

    #################################
    # simulation & save fig
    #################################
    time, X = maneuver_u_fix(K, T, x0, y0, psi0, u0, r0, duration, sampling, delta_list)
    draw_trajectory([X.T[1], X.T[2]], 'trajectory.png')
    draw_parameter_graph(time, X.T[0], X.T[3], X.T[4], X.T[5], 'parameter.png')
