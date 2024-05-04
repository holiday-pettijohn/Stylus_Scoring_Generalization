# Installing Stylus Workflow on the Compsci-ai Server

## Configuring Server Access

You will need to request access to the server manually, unless your research faculty have already configured your account. Contact the IT department at it.helpdesk@biola.edu to request this. In either case, you will have a personal username and password which you will need to access the server.

## Configuring VPN Access

To gain access to Biola's VPN which will in turn allow you to access the `compsci-ai` server from remote/home networks, ask Dr. Buzi to have your student account registered with Biola's campus VPN service. Biola VPN uses Cisco AnyConnect. Once you recieve information about the university VPN server, you will be able to use their tools to connect to the VPN and access the server from anywhere.

## Configuring Python Dependencies

## Installing Stylus

The libraries neccesary to interface with Stylus directly are located [here](https://github.com/biologic/stylus).

Open a shell on the server and run the following:

```
source /opt/anaconda/bin/activate root
git clone https://github.com/biologic/stylus.git
cd stylus
python setup.py install
```

You may need to run `python setup.py install` as sudo

This will activate your Anaconda enviornment and install the Python library for Stylus for your user

If the installation was successful, you should be able to run `python stylus_test.py` to ensure you can load the Stylus engine succesfully. If this command returns without error, you should be able to use the rest of the Stylus-related libraries included in this repository.

## Using the Installed Libraries

Whenever you plan to work with the Stylus libraries, be sure to activate your Anaconda enviorment with the following line:

`source /opt/anaconda/bin/activate root`
