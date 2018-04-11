# ship_mmg
"ship_mmg" is a web-based simulator based on Maneuvering Modeling Group (MMG) model for ship maneuvering.

## Description
This software can simulate the state of ship maneuvering from the time-series information of rudder angle(&delta;[rad]).

Example of the state of ship maneuvering is as follows:
- &Psi; [rad]
- Velocity and acceleration by ship coordinate system (X[m], Y[m], r[rad/s], etc.)
- etc.

You can change the simulation condition and the specification of target ship. In addition, you can get the simulation result of ship maneuvering by CSV data from web.
This is a web-based application by using Flask and Python. Simulator is implemented based on Python, Numpy and Scipy.

## Demo
![Interface](https://github.com/taiga4112/ship_mmg/wiki/images/demo_readme.png "Interface")

## Usage
Please asscess index page.

## Install and test on web browser
#### Mac or Linux
1. Fork it ([http://github.com/mitsuyukiLab/ship_mmg/fork](http://github.com/taiga4112/ship_mmg/fork))

2. Set developing environment
	```bash
	$ cd 'yourworkspace'
	$ git clone git@github.com:youraccount/ship_mmg.git
	$ virtualenv ship_mmg
	$ source ship_mmg/bin/activate
	$ pip install Flask numpy scipy matplotlib
	```

3. Start Flask app
	```bash
	$ python ship_mmg/view/app.py
	```

4. [Access](http://localhost:5000/)


#### Windows (installed by Anaconda)
1. Fork it ([http://github.com/mitsuyukiLab/ship_mmg/fork](http://github.com/mitsuyukiLab/ship_mmg/fork))

2. Set developing environment
  ```bash
  $ dir 'yourworkspace'
  $ git clone git@github.com:youraccount/ship_mmg.git
  $ conda create --name ship_mmg
  $ activate ship_mmg
  $ conda install Flask numpy scipy pandas matplotlib
  ```

3. Start Flask app
  ```bash
  $ python ship_mmg/view/app.py
  ```
4. [Access](http://localhost:5000/)


### for Developer and Researcher
Our laboratory has a research topic of creating digital twin of ship by using IoT data during voyage and system big data analysis technology based on this repository. Please contact us if you are interested in this topic.

mitsuyuki-taiga-my  at-mark  ynu.ac.jp


## Contribution
1. Fork it ( http://github.com/mitsuyukiLab/ship_mmg/fork )
2. Create your feature branch (git checkout -b my-new-feature)
3. Commit your changes (git commit -am 'Add some feature')
4. Push to the branch (git push origin my-new-feature)
5. Create new Pull Request

## Author

[taiga4112](https://github.com/taiga4112)
