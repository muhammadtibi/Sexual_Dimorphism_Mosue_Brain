# matlab_settings_gui
A small framework to simply ask input from users in a GUI. It can be considered as an extension to the built-in inputdlg of Matlab as it supports not only text input but many other types as well.

# Motivation
This small toolbox was designed to simplify parameter inqueries from users. It aims to eliminate the need of creating separate GUIs with customized layouts each and every time one needs to ask for parameters from users. Instead, by a unified parameter description interface it automatically generates the required GUI.

# Usage

The simplest usage of the toolbox is via the Settings_GUI function (GUIDE figure). You only need to create the appropriate parameter descirption variable (paramarray, see description in the doc of generateUIControls.m) and a modal figure is automatically created (i.e. it blocks the running of other code until the user enters the information). Once the user filled in the values and click OK the gui returns with a cellarray containing the specified value (for the exact format please check the documentation of fetchUIControlValues.m)

It is possible to provide a check function handle to the GUI so that it only accepts values that pass through the evaluation carried out by this check function.

# Example

Define a parameter structure and then call the Settings_GUI file

```
S = {...
  struct('name','An enum','type','enum','values',{{'Option #1','Option #2','Option #3'}});...
  struct('name','A number','type','int','default',5);...
  struct('name','A string','type','str','default','Sample text');...
  struct('name','Your favourite color','type','colorPicker','default',[0 1 0.8]);...
};

answers = Settings_GUI(S)
```
