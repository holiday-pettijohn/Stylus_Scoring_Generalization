# Running Jupyter on Compsci-ai Server

## Running Jupyter Remotely

In order to interface conveniently with the server enviornment you have configured, you may use a remote Jupyter notebook which is connected directly to the server. You will need to be on the Biola network in order to access the server. Refer to the section on 'Configuring VPN Access' if you are accessing the server from off campus.

This section assumes you have previously configured server access. If you have not, refer to 'Configuring Server Access'.

SSH into the compsci-ai server:

`ssh USERNAME@compsci-ai.biola.edu`

You will be prompted for your password, after which you will be presented with a shell on the server. Change your working directory to the one of your preference. Given the default configurations of this project, it is best to run Jupyter from the root directory of the project repository cloned from Git, wherever you have cloned it. For instance, `cd ~/Stylus_Scoring_Generalization`, if you cloned the code into your home directory.

Run the following command on the server in order to open a Jupyter notebook instance. Communicate with other users to ensure you choose a port which is not currently in use. Port `8890` is often used as a default, though assigning a different port for each person who is likely to need to use the notebook in paralell is reccomended.

`jupyter notebook --port YOUR_PORT_HERE --allow-root --ip=*`

Refer to the terminal output on this command and search for a link with the following structure:

`http://compsci-ai:YOUR_PORT_HERE/token=RANDOM_TOKEN_HERE`

You will need to edit this link by replacing 'compsci-ai' with 'compsci-ai.biola.edu' as below:

`http://compsci-ai.biola.edu:YOUR_PORT_HERE/token=RANDOM_TOKEN_HERE`

You should then be able to connect to the notebook via that link. Be aware that your notebook will start in whatever directory you originally ran the `jupyter notebook` command. You can run the command in any directory you choose, but you will not be able to easily access the parent folder of the directory you decide to run in.
