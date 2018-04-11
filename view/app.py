#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import sys
# reload(sys)
# sys.setdefaultencoding("utf-8")

import sys,os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
print(sys.path)
import kt_maneuver

from flask import Flask, render_template, request, redirect, url_for
import numpy as np

project_name = 'Ship MMG'
error_message = 'Value error. Please try again.'

app = Flask(__name__) # Flask setting


#############################################
# Routing setting for web application
#############################################

@app.route('/')
def index():
    title = project_name
    return render_template('index.html', title=title) # rendering 'index.html'


@app.route('/kt_fix', methods=['Get', 'Post'])
def kt_fix_view():
    title = project_name + " (KT View)"

    if request.method == 'POST':
        K = request.form.get('K', 0.155)
        T = request.form.get('T', 80.5)
        u0 = request.form.get('u0', 10.28)  # [m/s]
        x0 = request.form.get('x0', 0.0)  # [m]
        y0 = request.form.get('y0', 0.0)  # [m]
        psi0_degree = request.form.get('psi0_degree', 0.0)  # [degree]
        r0_degree = request.form.get('r0_degree', 0.0)  # [degree/s]
        duration = request.form.get('duration', 500)  # [s]
        sampling = request.args.get('sampling', 1000)  # [-]
        delta_degree = request.form.get('delta_degree', 5.0)  # [degree]

        shipL = request.form.get('shipL', 394.0)
        shipB = request.form.get('shipB', 50.6)
        shipD = request.form.get('shipD', 29.5)
        animation_speed = request.form.get('animation_speed', 10)

        # Unit conversion
        psi0 = float(psi0_degree) * np.pi / 180
        r0 = float(r0_degree) * np.pi / 180

        # KT simulation
        delta_list = np.ones(int(sampling)) * float(delta_degree) * np.pi / 180  # [rad]
        time, X = kt_maneuver.maneuver_u_fix(float(K), float(T), float(x0), float(y0), float(psi0), float(u0),
                                             float(r0), float(duration), float(sampling), delta_list)
        # Error check has to be done. But not yet implemented...

        return render_template('kt_fix_view.html', title=title, result=zip(time, X), K=K, T=T, u0=u0, x0=x0, y0=y0,
                               psi0_degree=psi0_degree, r0_degree=r0_degree, duration=duration, sampling=sampling,
                               delta_degree=delta_degree, shipL=shipL, shipD=shipD, shipB=shipB,
                               animation_speed=animation_speed)

    else:
        K = request.args.get('K', 0.155)
        T = request.args.get('T', 80.5)
        u0 = request.args.get('u0', 10.28)  # [m/s]
        x0 = request.args.get('x0', 0.0)  # [m]
        y0 = request.args.get('y0', 0.0)  # [m]
        psi0_degree = request.args.get('psi0_degree', 0.0)  # [degree]
        r0_degree = request.args.get('r0_degree', 0.0)  # [degree/s]
        duration = request.args.get('duration', 500)  # [s]
        sampling = request.args.get('sampling', 1000)  # [-]
        delta_degree = request.args.get('delta_degree', 5.0)  # [degree]

        shipL = request.args.get('shipL', 364.0)
        shipB = request.args.get('shipB', 50.6)
        shipD = request.args.get('shipD', 29.5)
        animation_speed = request.args.get('animation_speed', 10)

        return render_template('kt_fix_view.html', title=title, K=K, T=T, u0=u0, x0=x0, y0=y0, psi0_degree=psi0_degree,
                               r0_degree=r0_degree, duration=duration, sampling=sampling, delta_degree=delta_degree,
                               shipL=shipL, shipD=shipD, shipB=shipB, animation_speed=animation_speed)


@app.route('/kt_trigomestric', methods=['Get', 'Post'])
def kt_trigomestric_view():
    title = project_name + " (KT View)"

    if request.method == 'POST':
        K = request.form.get('K', 0.155)
        T = request.form.get('T', 80.5)
        u0 = request.form.get('u0', 10.28)  # [m/s]
        x0 = request.form.get('x0', 0.0)  # [m]
        y0 = request.form.get('y0', 0.0)  # [m]
        psi0_degree = request.form.get('psi0_degree', 0.0)  # [degree]
        r0_degree = request.form.get('r0_degree', 0.0)  # [degree/s]
        duration = request.form.get('duration', 500)  # [s]
        sampling = request.form.get('sampling', 1000)  # [-]
        delta0_degree = request.form.get('delta0_degree', 10.0)  # [degree]
        w0_degree = request.form.get('w0_degree', 14.4)  # [degree/s]
        theta0_degree = request.form.get('theta0_degree', 0.0)  # [degree]

        shipL = request.form.get('shipL', 394.0)
        shipB = request.form.get('shipB', 50.6)
        shipD = request.form.get('shipD', 29.5)
        animation_speed = request.form.get('animation_speed', 10)

        # Unit conversion
        psi0 = float(psi0_degree) * np.pi / 180
        r0 = float(r0_degree) * np.pi / 180

        # KT simulation
        time_list = np.linspace(0.00, float(duration), float(sampling))
        _degree = float(w0_degree) * np.pi / 180.0 * time_list + float(theta0_degree) * np.pi / 180.0
        delta_list = np.sin(_degree) * float(delta0_degree) * np.pi / 180  # [rad]
        time, X = kt_maneuver.maneuver_u_fix(float(K), float(T), float(x0), float(y0), float(psi0), float(u0),
                                             float(r0), float(duration), float(sampling), delta_list)
        # Error check has to be done. But not yet implemented...

        return render_template('kt_trigomestric_view.html', title=title, result=zip(time, X), K=K, T=T, u0=u0, x0=x0,
                               y0=y0, psi0_degree=psi0_degree, r0_degree=r0_degree, duration=duration,
                               sampling=sampling, delta0_degree=delta0_degree, w0_degree=w0_degree,
                               theta0_degree=theta0_degree, shipL=shipL, shipD=shipD, shipB=shipB,
                               animation_speed=animation_speed)

    else:
        K = request.args.get('K', 0.155)
        T = request.args.get('T', 80.5)
        u0 = request.args.get('u0', 10.28)  # [m/s]
        x0 = request.args.get('x0', 0.0)  # [m]
        y0 = request.args.get('y0', 0.0)  # [m]
        psi0_degree = request.args.get('psi0_degree', 0.0)  # [degree]
        r0_degree = request.args.get('r0_degree', 0.0)  # [degree/s]
        duration = request.args.get('duration', 500)  # [s]
        sampling = request.args.get('sampling', 1000)  # [-]
        delta0_degree = request.args.get('delta0_degree', 10.0)  # [degree]
        w0_degree = request.args.get('w0_degree', 14.4)  # [degree/s]
        theta0_degree = request.args.get('theta0_degree', 0.0)  # [degree]

        shipL = request.args.get('shipL', 364.0)
        shipB = request.args.get('shipB', 50.6)
        shipD = request.args.get('shipD', 29.5)
        animation_speed = request.args.get('animation_speed', 10)

        return render_template('kt_trigomestric_view.html', title=title, K=K, T=T, u0=u0, x0=x0, y0=y0,
                               psi0_degree=psi0_degree, r0_degree=r0_degree, duration=duration, sampling=sampling,
                               delta0_degree=delta0_degree, w0_degree=w0_degree, theta0_degree=theta0_degree,
                               shipL=shipL, shipD=shipD, shipB=shipB, animation_speed=animation_speed)


@app.route('/kt_zigzag', methods=['Get', 'Post'])
def kt_zigzag_view():
    title = project_name + " (KT View)"

    if request.method == 'POST':
        K = request.form.get('K', 0.155)
        T = request.form.get('T', 80.5)
        u0 = request.form.get('u0', 10.28)  # [m/s]
        x0 = request.form.get('x0', 0.0)  # [m]
        y0 = request.form.get('y0', 0.0)  # [m]
        psi0_degree = request.form.get('psi0_degree', 0.0)  # [degree]
        r0_degree = request.form.get('r0_degree', 0.0)  # [degree/s]
        duration = request.form.get('duration', 500)  # [s]
        sampling = request.form.get('sampling', 1000)  # [-]
        delta0_degree = request.form.get('delta0_degree', -10.0)  # [degree]
        psi_a_degree = request.form.get('psi_a_degree', 20.0)  # [degree/s]
        t1 = request.form.get('t1', 1.0)  # [degree]

        shipL = request.form.get('shipL', 394.0)
        shipB = request.form.get('shipB', 50.6)
        shipD = request.form.get('shipD', 29.5)
        animation_speed = request.form.get('animation_speed', 10)

        # Unit conversion
        psi0 = float(psi0_degree) * np.pi / 180
        r0 = float(r0_degree) * np.pi / 180

        # KT simulation
        delta_a = float(delta0_degree) * np.pi / 180.0
        psi_a = float(psi_a_degree) * np.pi / 180.0
        time, X = kt_maneuver.maneuver_zigzag(float(K), float(T), float(x0), float(y0), float(psi0), float(u0),
                                              float(r0), float(duration), float(sampling), delta_a, psi_a, float(t1))
        # Error check has to be done. But not yet implemented...

        return render_template('kt_zigzag_view.html', title=title, result=zip(time, X), K=K, T=T, u0=u0, x0=x0, y0=y0,
                               psi0_degree=psi0_degree, r0_degree=r0_degree, duration=duration, sampling=sampling,
                               delta0_degree=delta0_degree, psi_a_degree=psi_a_degree, t1=t1, shipL=shipL, shipD=shipD,
                               shipB=shipB, animation_speed=animation_speed)

    else:
        K = request.args.get('K', 0.155)
        T = request.args.get('T', 80.5)
        u0 = request.args.get('u0', 10.28)  # [m/s]
        x0 = request.args.get('x0', 0.0)  # [m]
        y0 = request.args.get('y0', 0.0)  # [m]
        psi0_degree = request.args.get('psi0_degree', 0.0)  # [degree]
        r0_degree = request.args.get('r0_degree', 0.0)  # [degree/s]
        duration = request.args.get('duration', 500)  # [s]
        sampling = request.args.get('sampling', 1000)  # [-]
        delta0_degree = request.args.get('delta0_degree', -10.0)  # [degree]
        psi_a_degree = request.args.get('psi_a_degree', 20.0)  # [degree/s]
        t1 = request.args.get('t1', 1.0)  # [degree]

        shipL = request.args.get('shipL', 364.0)
        shipB = request.args.get('shipB', 50.6)
        shipD = request.args.get('shipD', 29.5)
        animation_speed = request.args.get('animation_speed', 10)

        return render_template('kt_zigzag_view.html', title=title, K=K, T=T, u0=u0, x0=x0, y0=y0,
                               psi0_degree=psi0_degree, r0_degree=r0_degree, duration=duration, sampling=sampling,
                               delta0_degree=delta0_degree, psi_a_degree=psi_a_degree, t1=t1, shipL=shipL, shipD=shipD,
                               shipB=shipB, animation_speed=animation_speed)


# Get angle(rad) from angle(degree)
def __rad(delta_degree):
    delta_rad = float(delta_degree) * np.pi / 180
    return delta_rad


#######################################
# main setting
#######################################
if __name__ == '__main__':
    # app.debug = True # for DEBUG
    app.run(host='0.0.0.0') # for all access
