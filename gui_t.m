function varargout = gui_t(varargin)
% GUI_T MATLAB code for gui_t.fig
%      GUI_T, by itself, creates a new GUI_T or raises the existing
%      singleton*.
%
%      H = GUI_T returns the handle to a new GUI_T or the handle to
%      the existing singleton*.
%
%      GUI_T('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_T.M with the given input arguments.
%
%      GUI_T('Property','Value',...) creates a new GUI_T or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_t_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_t_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_t

% Last Modified by GUIDE v2.5 03-Sep-2021 18:33:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_t_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_t_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_t is made visible.
function gui_t_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_t (see VARARGIN)

% Choose default command line output for gui_t
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_t wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_t_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in start_but.
function start_but_Callback(hObject, eventdata, handles)
% hObject    handle to start_but (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(hObject,'Enable','off');
set(handles.clear_but,'Enable','on');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.func_choice,'Enable','on');
set(handles.cp_fun,'Enable','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Extracting data from the table  %
%  with the physical parameters    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = get(handles.table_val,'Data');
sz = size(data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m               =    str2num(data{1,1});
handles.m       =    m;

m0              =    str2num(data{2,1});
handles.m0      =    m0;

cs              =    str2num(data{3,1});
handles.cs      =    cs;

ms              =    str2num(data{4,1});
handles.ms      =    ms;

S               =    str2num(data{5,1});
handles.S       =    S;

sigma           =    str2num(data{6,1});
handles.sigma   =    sigma;

k               =    str2num(data{7,1});
handles.k       =    k;

eps             =    str2num(data{8,1});
handles.eps     =    eps;

t0              =    str2num(data{9,1});
handles.t0      =    t0;

tend            =    str2num(data{10,1});
handles.tend    =    tend;

yinit           =    str2num(data{11,1});
handles.yinit   =    yinit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Saving handles structure     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Extracting data from the table  %
%        with t, T(t), M(t)        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dat = get(handles.table_data,'Data');
sz = size(Dat);

for i = 1 : sz(1)
    t(i) = Dat(i,1);
    T(i) = Dat(i,2);
    M(i) = Dat(i,3);
end

handles.t = t; % dimensional t!!!
handles.T = T;
handles.M = M;
tab_dat = get(handles.table_data,'Data');
sz_tab_dat = size(tab_dat);
handles.num_data = ...
    sz_tab_dat(1);%numel(handles.t);
%display(handles.num_data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Saving handles structure     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracting start and final index %
%             for time             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start_s         =       get(handles.start_index,'String'); 
%start_i         =       str2num(start_s);

%end_s           =       get(handles.end_index,'String'); 
%end_i           =       str2num(end_s);

handles.v = handles.start_i : handles.end_i;
%handles.v = v;

%%%%%%%%%%%%%%%%%%%%%%
% DIMENSIONLESS GRID %
%%%%%%%%%%%%%%%%%%%%%%
t0_dimless = handles.t(handles.start_i)/handles.tend;%handles.t(handles.start_i);
handles.t0_dimless = t0_dimless;
tend_dimless = 1;
handles.tend_dimless = tend_dimless;
t_dimless = ...
    linspace(t0_dimless,tend_dimless,numel(handles.t(handles.start_i : handles.end_i)));
handles.t_dimless = t_dimless;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Saving handles structure     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i = 1 : sz(1)
%    col(i) = str2num(data{i,1});
%end
%s = 0;
%for i = 1 : length(col)
%    s = s + col(i);
%end
%f = warndlg(['Sum: ',num2str(s)],'Information');
%f = msgbox(['Size: ',num2str(size(data))],'User information','help');
%f = warndlg(['Element: ',num2str(data{2,1})],'Information');
%f = warndlg(['col: ',num2str(col)],'Information');
%f = warndlg(['mean(col): ',num2str(mean(col))],'Information');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dat = load('data_1.txt');
%t = Dat(:,1); 
%T = Dat(:,2); 
%M = Dat(:,3); 
%d = get(handles.table_data,'Data');
%for i = 1 : length(t)
%    d{i,1} = t(i);
%    d{i,2} = T(i);
%    d{i,3} = M(i);
%end
%p = warndlg(['Size of table_data: ',num2str(size(d))],'Information');
%t0_ind_s = get(handles.start_index,'String'); t0_ind = str2num(t0_ind_s);
%tend_ind_s = get(handles.end_index,'String'); tend_ind = str2num(tend_ind_s);

% --- Executes on button press in clear_but.
function clear_but_Callback(hObject, eventdata, handles)
% hObject    handle to clear_but (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tab_dat = get(handles.table_data,'Data');
sz_tab_dat = size(tab_dat);
ndata = sz_tab_dat(1);

set(hObject,'Enable','off');

set(handles.start_index,'Enable','on');
set(handles.start_index,'String','1');

set(handles.end_index,'Enable','on');
%set(handles.end_index,'String','326');
set(handles.end_index,'String',num2str(ndata))

set(handles.func_choice,'Enable','off');
set(handles.cp_fun,'Enable','off');

set(handles.cp_fun,'Value',1);
set(handles.func_choice,'Value',1);

axes(handles.ax_M_form)
set(handles.ax_M_form,'Visible','off')
%delete(handles.MText);
cla(handles.ax_M_form);

axes(handles.ax_c0_form)
set(handles.ax_c0_form,'Visible','off');
%delete(handles.c0Text);
cla(handles.ax_c0_form);

sz = size(get(handles.table_val,'Data'));
vec = {'';'';'';'';'';'';'';'';'';'';''};
set(handles.table_val,'Data',vec);
%class(get(handles.table_val,'Data'))

axes(handles.ax_main)
%set(handles.hLeg_data_plot,'Visible','off');
set(handles.hLeg_num_plot,'Visible','off');
grid off;
cla;
%handles = rmfield(handles,'data_plot');
%handles = rmfield(handles,'num_plot');
%handles = rmfield(handles,'hLeg_data_plot');
%handles = rmfield(handles,'hLeg_num_plot');

axes(handles.c_ax)
set(handles.hLeg_cp,'Visible','off');
grid off;
cla;
%handles = rmfield(handles,'cplot');
%handles = rmfield(handles,'hLeg_cp');

set(handles.est_vals,'Visible','off')
        
set(handles.p1_val,'Visible','off')
set(handles.p2_val,'Visible','off')
set(handles.p3_val,'Visible','off')
set(handles.p4_val,'Visible','off')
set(handles.p5_val,'Visible','off')
        
set(handles.p1_ed,'Visible','off')
set(handles.p1_ed,'String','')
set(handles.p2_ed,'Visible','off')
set(handles.p2_ed,'String','')
set(handles.p3_ed,'Visible','off')
set(handles.p3_ed,'String','')
set(handles.p4_ed,'Visible','off')
set(handles.p4_ed,'String','')
set(handles.p5_ed,'Visible','off')
set(handles.p5_ed,'String','')

set(handles.export_temp,'Enable','off')
set(handles.export_heat,'Enable','off')
set(handles.plot_heat,'Enable','off')

        

% --- Executes when entered data in editable cell(s) in table_val.
function table_val_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_val (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)



% --- Executes when selected cell(s) is changed in table_val.
function table_val_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_val (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


%{
data = get(handles.table_val,'Data');
len = length(data);
for i = 1 : len
    d(i) = str2num(data{i});
end
res = true;
for i = 1 : length(d)
    if(~isnumeric(d(i)))
        res = false;
    end
end
% START button activation
if(res==true)
    set(handles.start_but,'Enable','on');
end
%}


% --- Executes on key press with focus on table_val and none of its controls.
function table_val_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to table_val (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function table_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to table_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --------------------------------------------------------------------
function table_val_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to table_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extracting physical parameters values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data = get(hObject,'Data');
data = get(handles.table_val,'Data');
sz = size(data);
res = true;
handles.num_points = sz(1);
guidata(hObject,handles);
%class(data)

for i = 1 : sz(1)
    if(isempty(data{i,1}) || ...
       ~isnumeric(str2num(data{i,1})) || ...
       strcmp(num2str(str2num(data{i,1})),''))
        res = false;
    end
end
%}
%{
for i = 1 : sz(1)
    if(isempty(data(i,1)) || ...
       ~isnumeric(data(i,1)) || ...
       strcmp(num2str(data(i,1)),''))
        res = false;
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START button activation %
% at certain condition    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode_start = get(handles.start_index,'Enable');
res_start = strcmp(mode_start,'off');

mode_end = get(handles.end_index,'Enable');
res_end = strcmp(mode_end,'off');

if(res==true && res_start==true && res_end==true)
    set(handles.start_but,'Enable','on');
end


% --- Executes on selection change in func_choice.
function func_choice_Callback(hObject, eventdata, handles)
% hObject    handle to func_choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns func_choice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from func_choice

val = get(hObject,'Value');
switch val
    case 1
        %=============================================
        % c0
        %=============================================
        set(handles.ax_c0_form,'Visible','on');
        axes(handles.ax_c0_form)
        pos = get(handles.ax_c0_form,'position');
        handles.c0Text = ...
            text(pos(3)/2,pos(4)/2,'$$ c_{0}(T) = e\frac{aT+b}{cT+d} $$', ...
            'interpreter','latex',...
            'FontSize',9);
        set(handles.ax_c0_form,'TickDir','out');
        set(handles.ax_c0_form,'visible','off');
        axis([0,pos(3),0,pos(4)])
        %=============================================
        % input dialogue box
        %=============================================
        x = inputdlg({'a:','b:','c:','d:','e:'}, ...
            'Coefficients of c0(T)', ...
            [1,75; 1,75; 1,75; 1,75; 1,75]);
        %=============================================
        % Extract input values for c0(T)
        %=============================================
        handles.c0_a = str2num(x{1});
        handles.c0_b = str2num(x{2});
        handles.c0_c = str2num(x{3});
        handles.c0_d = str2num(x{4});
        handles.c0_e = str2num(x{5});
        handles.c0_args = ...
            [handles.c0_a, handles.c0_b, ...
             handles.c0_c, handles.c0_d, ...
             handles.c0_e];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Saving handles structure     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        guidata(hObject,handles);
        
    case 2
        %=============================================
        % M
        %=============================================
        set(handles.ax_M_form,'Visible','on')
        axes(handles.ax_M_form)
        pos = get(handles.ax_M_form,'position');
        handles.MText = ... 
            text(pos(3)/2,pos(4)/2,'$$ M(t) = ae^{bt} + ce^{dt} $$', ...
            'interpreter','latex',...
            'FontSize',9);
        set(handles.ax_M_form,'TickDir','out');
        set(handles.ax_M_form,'visible','off');
        axis([0,pos(3),0,pos(4)])
        %=============================================
        % input dialogue box
        %=============================================
        x = inputdlg({'a:','b:','c:','d:'}, ...
            'Coefficients of M(t)', ...
            [1,75; 1,75; 1,75; 1,75]);
        %=============================================
        % Extract input values for M(t)
        %=============================================
        handles.M_a = str2num(x{1});
        handles.M_b = str2num(x{2});
        handles.M_c = str2num(x{3});
        handles.M_d = str2num(x{4});
        handles.M_args = ...
            [handles.M_a, handles.M_b, ...
             handles.M_c, handles.M_d];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Saving handles structure     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        guidata(hObject,handles);
end


% --- Executes during object creation, after setting all properties.
function func_choice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to func_choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=================================
% M(t)
%=================================
function res = M_fun(M_args,time)
    a = M_args(1);
    b = M_args(2);
    c = M_args(3);
    d = M_args(4);
    res = a * exp(b * time) + ...
          c * exp(d*time);
%end

%=================================
% c0(T)
%=================================
function res = c0_fun(c0_args,T)
    a = c0_args(1);
    b = c0_args(2);
    c = c0_args(3);
    d = c0_args(4);
    e = c0_args(5);
    s = -1;
    numer = a * (T + s * 273.15) + b;
    denom = c * (T + s * 273.15) + d;
    res = e * (numer ./ denom);
%end

%=================================
% Polynomial c(T) 
%=================================
function res = cp_eval(p,u)
    n = length(p);
    res = 0;
    for i = 1 : n
        res = res + p(i)*u.^(i-1);
    end
%end

% --- Executes on selection change in cp_fun.
function cp_fun_Callback(hObject, eventdata, handles)
% hObject    handle to cp_fun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cp_fun contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cp_fun

axes(handles.ax_c0_form)
cla(handles.ax_c0_form);

axes(handles.ax_M_form)
cla(handles.ax_M_form);

val = get(hObject,'Value');
switch val
    case 1
        axes(handles.ax_main)
        cla(handles.ax_main)
        %set(handles.hLeg_data_plot,'Visible','off');
        %set(handles.hLeg_num_plot,'Visible','off');
        
        axes(handles.c_ax)
        cla(handles.c_ax)
        %set(handles.hLeg_cp,'Visible','off');
        %=============================================
        % input dialogue box, 1-st degree c(T)
        %=============================================
        %set(hObject,'Enable','off');
        %%%%%%%%%%%%%%%%
        %    Bounds    %
        %%%%%%%%%%%%%%%%
        bnd_prompt = {'Lower bounds for [p_1,p_2]', ...
          'Upper bounds for [p_1,p_2]'};
        bnd_wtitle = 'Coefficient bounds for c(T) = p_{2}*T + p_{1}';
        bnd_dims = [1,85];
        bnd_definput = {'',''};
        opts.Interpreter = 'tex';
        bnd = ...
            inputdlg(bnd_prompt,bnd_wtitle,bnd_dims,bnd_definput,opts);
        lb = str2num(bnd{1});
        handles.lb = lb;
        ub = str2num(bnd{2});
        handles.ub = ub;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check for correct input values %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while(numel(lb)~=2 || numel(ub)~=2)
            err = errordlg(['Out of bounds. Please try with correct ', ...
                'number of input values.'],'Input error');
            pause(3);
            %close(err);
            bnd_prompt = {'Lower bounds for [p_1,p_2]', ...
                  'Upper bounds for [p_1,p_2]'};
            bnd_wtitle = 'Coefficient bounds for c(T) = p_{2}*T + p_{1}';
            bnd_dims = [1,85];
            bnd_definput = {'',''};
            opts.Interpreter = 'tex';
            bnd = inputdlg(bnd_prompt,bnd_wtitle,bnd_dims,bnd_definput,opts);
            lb = str2num(bnd{1});
            ub = str2num(bnd{2});
        end
        %%%%%%%%%%%%%%%%%%%%%
        % Extracting bounds %
        %%%%%%%%%%%%%%%%%%%%%
        lb = str2num(bnd{1});
        handles.lb = lb;
        ub = str2num(bnd{2});
        handles.ub = ub;
        %%%%%%%%%%%%%%%%
        % Coefficients %
        %%%%%%%%%%%%%%%%
        prompt = {'p_{1}','p_{2}'};
        wtitle = ...
            'Initial coefficients of c(T) = p_{2}*T + p_{1}';
        dims = [1,85];
        definput = {'',''};
        opts.Interpreter = 'tex';
        x = inputdlg(prompt,wtitle,dims,definput,opts);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Initial values for p0  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        p0(1) = str2num(x{1});
        p0(2) = str2num(x{2});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Estimating p vector    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %handles.fun = @(t,u) 2*u+3*t;
        %w_info = msgbox('Function f(t,u) successfully created.', ...
        %                'User information','help');
        %p = p0; %[1,1];
        %ydata = handles.T;
        anon_f_1 = @(p,tdata) fun_deg_1(p,tdata,handles);
        tdata = handles.t_dimless; %linspace(0,1,numel(handles.t));
        ydata = handles.T(handles.start_i:handles.end_i); %fun_deg_1(p0,tdata,handles);
        handles.p_est = lsqcurvefit(anon_f_1,p0,tdata,ydata,lb,ub);
        p_win = msgbox('Parameter values for 1-st degree c(T) estimated.', ...
            'User information','help');
        
        Max_err = max(abs(ydata-fun_deg_1(handles.p_est,tdata,handles)));
        Max_err_win = msgbox(['Maximal absoulte error |y-u|: ', ...
            num2str(Max_err)],'User information','help');
        
        set(handles.est_vals,'Visible','on')
        
        set(handles.p1_val,'Visible','on')
        set(handles.p2_val,'Visible','on')
        
        set(handles.p1_ed,'Visible','on')
        set(handles.p1_ed,'String',handles.p_est(1))
        set(handles.p2_ed,'Visible','on')
        set(handles.p2_ed,'String',handles.p_est(2))
        
        set(handles.export_temp,'Enable','on')
        set(handles.export_heat,'Enable','on')
        set(handles.plot_heat,'Enable','on')
        %r = fun_deg_1(p0,tdata,handles);
        %win = msgbox('The handles structure passed as an argument.', ...
        %                'User information','help');
        %r_win = msgbox('Numerical solution with 1-st degree c(T) found.', ...
        %                'User information','help');
        %fields(handles)
        %%%%%%%%%%%%%%%%%%%%%%%
        % Plot estimated c(T) %
        %%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.c_ax)
        handles.cplot = ...
            plot(handles.T-273.15,cp_eval(handles.p_est,handles.T));
        set(handles.cplot,'LineWidth',3,'LineStyle',':','Color','r')
        set(handles.c_ax,'FontSize',8)
        grid on
        xlabel('\bf{T}')
        ylabel('\bf{c(T)}')
        handles.hLeg_cp = legend('\it{1-st degree c(T) approximation}');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Connecting with the context menu %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.cplot,'UIContextMenu',handles.line_properties)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot dimensional solution with % 
        % estimated parameters and the   %
        % experimental data  T(t)        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.ax_main)
        %cla(handles.ax_main)
        set(handles.ax_main,'FontSize',8)
        hold on
        grid on
        xlabel('\bf{t}')
        ylabel('\bf{T(t), [C deg]}')
        %%%%%%%%%%%%%%%%%%%%%
        % Experimental data %
        %%%%%%%%%%%%%%%%%%%%%
        handles.data_plot = ...
            plot(handles.t(handles.v), ...
            handles.T(handles.v)-273.15);
        set(handles.data_plot,'Color','b', ...
                              'LineStyle','-', ...
                              'LineWidth',3)
        %handles.hLeg_data_plot = ...
        %    legend('\it{Experimental data}');
        %%%%%%%%%%%%%%%%%%%%%%
        % Numerical solution %
        %%%%%%%%%%%%%%%%%%%%%%
        handles.num_sol = ...
            fun_deg_1(handles.p_est,tdata,handles);
        %axes(handles.ax_main)
        %hold on
        handles.num_plot = ...
            plot(handles.t,handles.num_sol-273.15);
        set(handles.num_plot,'Color','g', ...
                             'LineStyle',':', ...
                             'LineWidth',3)
        handles.hLeg_num_plot = ...
            legend('\it{Experimental data}','\it{Numerical solution}');
        %%%%%%%%%%%%%%%%%%%%%%
        %    Save handles    %
        %%%%%%%%%%%%%%%%%%%%%%
        guidata(hObject,handles);
    case 2
        axes(handles.ax_main)
        cla(handles.ax_main)
        
        axes(handles.c_ax)
        cla(handles.c_ax)
        %=============================================
        % input dialogue box, 2-nd degree c(T)
        %=============================================
        %set(hObject,'Enable','off');
        %%%%%%%%%%%%%%%%
        %    Bounds    %
        %%%%%%%%%%%%%%%%
        bnd_prompt = {'Lower bounds for [p_1,p_2,p_3]', ...
          'Upper bounds for [p_1,p_2,p_3]'};
        bnd_wtitle = 'Coefficient bounds for c(T) = p_{3}*T^2 + p_{2}*T + p_{1}';
        bnd_dims = [1,100];
        bnd_definput = {'',''};
        opts.Interpreter = 'tex';
        bnd = ...
            inputdlg(bnd_prompt,bnd_wtitle,bnd_dims,bnd_definput,opts);
        lb = str2num(bnd{1});
        handles.lb = lb;
        ub = str2num(bnd{2});
        handles.ub = ub;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check for correct input values %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while(numel(lb)~=3 || numel(ub)~=3)
            err = errordlg(['Out of bounds. Please try with correct ', ...
                'number of input values.'],'Input error');
            pause(3);
            %close(err);
            bnd_prompt = {'Lower bounds for [p_1,p_2,_3]', ...
                  'Upper bounds for [p_1,p_2,p_3]'};
            bnd_wtitle = ...
                'Coefficient bounds for c(T) = p_{3}*T^2 + p_{2}*T + p_{1}';
            bnd_dims = [1,100];
            bnd_definput = {'',''};
            opts.Interpreter = 'tex';
            bnd = inputdlg(bnd_prompt,bnd_wtitle,bnd_dims,bnd_definput,opts);
            lb = str2num(bnd{1});
            ub = str2num(bnd{2});
        end
        %%%%%%%%%%%%%%%%%%%%%
        % Extracting bounds %
        %%%%%%%%%%%%%%%%%%%%%
        lb = str2num(bnd{1});
        handles.lb = lb;
        ub = str2num(bnd{2});
        handles.ub = ub;
        %%%%%%%%%%%%%%%%
        % Coefficients %
        %%%%%%%%%%%%%%%%
        prompt = {'p_{1}','p_{2}','p_{3}'};
        wtitle = ...
            'Initial coefficients of c(T) = p_{3}*T^2 + p_{2}*T + p_{1}';
        dims = [1,95];
        definput = {'','',''};
        opts.Interpreter = 'tex';
        x = inputdlg(prompt,wtitle,dims,definput,opts);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Initial values for p0  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        p0(1) = str2num(x{1});
        p0(2) = str2num(x{2});
        p0(3) = str2num(x{3});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Estimating p vector    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %handles.fun = @(t,u) 2*u+3*t;
        %w_info = msgbox('Function f(t,u) successfully created.', ...
        %                'User information','help');
        %p = p0; %[1,1,1];
        %tdata = handles.t_dimless; %linspace(0,1,numel(handles.t));
        anon_f_2 = @(p,tdata) fun_deg_2(p,tdata,handles);
        %tdata = handles.t_dimless; %linspace(0,1,numel(handles.t));
        tdata = handles.t(handles.start_i:handles.end_i)/handles.t(end);
        ydata = handles.T(handles.start_i:handles.end_i); %fun_deg_2(p0,tdata,handles);
        handles.p_est = lsqcurvefit(anon_f_2,p0,tdata,ydata,lb,ub);
        p_win = msgbox('Parameter values for 2-nd degree c(T) estimated.', ...
            'User information','help');
        
        Max_err = max(abs(ydata-fun_deg_2(handles.p_est,tdata,handles)));
        Max_err_win = msgbox(['Maximal absoulte error |y-u|: ', ...
            num2str(Max_err)],'User information','help');
        
        set(handles.est_vals,'Visible','on')
        
        set(handles.p1_val,'Visible','on')
        set(handles.p2_val,'Visible','on')
        set(handles.p3_val,'Visible','on')
        
        set(handles.p1_ed,'Visible','on')
        set(handles.p1_ed,'String',handles.p_est(1))
        set(handles.p2_ed,'Visible','on')
        set(handles.p2_ed,'String',handles.p_est(2))
        set(handles.p3_ed,'Visible','on')
        set(handles.p3_ed,'String',handles.p_est(3))
        
        set(handles.export_temp,'Enable','on')
        set(handles.export_heat,'Enable','on')
        set(handles.plot_heat,'Enable','on')
        %r = fun_deg_2(p0,tdata,handles);
        %win = msgbox('The handles structure passed as an argument.', ...
        %                'User information','help');
        %r_win = msgbox('Numerical solution with 2-nd degree c(T) found.', ...
        %                'User information','help');
        %%%%%%%%%%%%%%%%%%%%%%%
        % Plot estimated c(T) %
        %%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.c_ax)
        handles.cplot = ...
            plot(handles.T-273.15,cp_eval(handles.p_est,handles.T));
        set(handles.cplot,'LineWidth',3,'LineStyle',':','Color','g')
        set(handles.c_ax,'FontSize',8)
        grid on
        xlabel('\bf{T}')
        ylabel('\bf{c(T)}')
        handles.hLeg_cp = legend('\it{2-nd degree c(T) approximation}');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot dimensional solution with % 
        % estimated parameters and the   %
        % experimental data  T(t)        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.ax_main)
        %cla(handles.ax_main)
        set(handles.ax_main,'FontSize',8)
        hold on
        grid on
        xlabel('\bf{t}')
        ylabel('\bf{T(t), [C deg]}')
        %%%%%%%%%%%%%%%%%%%%%
        % Experimental data %
        %%%%%%%%%%%%%%%%%%%%%
        handles.data_plot = ...
            plot(handles.t(handles.v), ...
            handles.T(handles.v)-273.15);
        set(handles.data_plot,'Color','b', ...
                              'LineStyle','-', ...
                              'LineWidth',3)
        handles.hLeg_data_plot = ...
            legend('\it{Experimental data}');
        %%%%%%%%%%%%%%%%%%%%%%
        % Numerical solution %
        %%%%%%%%%%%%%%%%%%%%%%
        handles.num_sol = ...
            fun_deg_2(handles.p_est,tdata,handles);
        %axes(handles.ax_main)
        %hold on
        handles.num_plot = ...
            plot(handles.t(handles.v),handles.num_sol-273.15);
        set(handles.num_plot,'Color','g', ...
                             'LineStyle',':', ...
                             'LineWidth',3)
        handles.hLeg_num_plot = ...
            legend('\it{Experimental data}','\it{Numerical solution}');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Connecting with the context menu %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.cplot,'UIContextMenu',handles.line_properties)
        %%%%%%%%%%%%%%%%%%%%%%
        %    Save handles    %
        %%%%%%%%%%%%%%%%%%%%%%
        guidata(hObject,handles);
    case 3
        axes(handles.ax_main)
        cla(handles.ax_main)
        
        axes(handles.c_ax)
        cla(handles.c_ax)
        %=============================================
        % input dialogue box, 3-rd degree c(T)
        %=============================================
        %set(hObject,'Enable','off');
        %%%%%%%%%%%%%%%%
        %    Bounds    %
        %%%%%%%%%%%%%%%%
        bnd_prompt = {'Lower bounds for [p_1,p_2,p_3,p_4]', ...
          'Upper bounds for [p_1,p_2,p_3,p_4]'};
        bnd_wtitle = 'Coefficient bounds for c(T) = p_{4}*T^3 + p_{3}*T^2 + p_{2}*T + p_{1}';
        bnd_dims = [1,115];
        bnd_definput = {'',''};
        opts.Interpreter = 'tex';
        bnd = ...
            inputdlg(bnd_prompt,bnd_wtitle,bnd_dims,bnd_definput,opts);
        lb = str2num(bnd{1});
        handles.lb = lb;
        ub = str2num(bnd{2});
        handles.ub = ub;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check for correct input values %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while(numel(lb)~=4 || numel(ub)~=4)
            err = errordlg(['Out of bounds. Please try with correct ', ...
                'number of input values.'],'Input error');
            pause(3);
            %close(err);
            bnd_prompt = {'Lower bounds for [p_1,p_2,_3,p_4]', ...
                  'Upper bounds for [p_1,p_2,p_3,p_4]'};
            bnd_wtitle = ...
                'Coefficient bounds for c(T) = p_{4}*T^3 + p_{3}*T^2 + p_{2}*T + p_{1}';
            bnd_dims = [1,115];
            bnd_definput = {'',''};
            opts.Interpreter = 'tex';
            bnd = inputdlg(bnd_prompt,bnd_wtitle,bnd_dims,bnd_definput,opts);
            lb = str2num(bnd{1});
            ub = str2num(bnd{2});
        end
        %%%%%%%%%%%%%%%%%%%%%
        % Extracting bounds %
        %%%%%%%%%%%%%%%%%%%%%
        lb = str2num(bnd{1});
        handles.lb = lb;
        ub = str2num(bnd{2});
        handles.ub = ub;
        %%%%%%%%%%%%%%%%
        % Coefficients %
        %%%%%%%%%%%%%%%%
        prompt = {'p_{1}','p_{2}','p_{3}','p_{4}'};
        wtitle = ...
            'Initial coefficients of c(T) = p_{4}*T^3 + p_{3}*T^2 + p_{2}*T + p_{1}';
        dims = [1,110];
        definput = {'','','',''};
        opts.Interpreter = 'tex';
        x = inputdlg(prompt,wtitle,dims,definput,opts);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Initial values for p0  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        p0(1) = str2num(x{1});
        p0(2) = str2num(x{2});
        p0(3) = str2num(x{3});
        p0(4) = str2num(x{4});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Estimating p vector    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %handles.fun = @(t,u) 2*u+3*t;
        %w_info = msgbox('Function f(t,u) successfully created.', ...
        %                'User information','help');
        %p = p0; %[1,1,1,1];
        %tdata = handles.t_dimless; %linspace(0,1,numel(handles.t));
        anon_f_3 = @(p,tdata) fun_deg_3(p,tdata,handles);
        tdata = handles.t_dimless; %linspace(0,1,numel(handles.t));
        ydata = handles.T(handles.start_i:handles.end_i); %fun_deg_3(p0,tdata,handles);
        handles.p_est = lsqcurvefit(anon_f_3,p0,tdata,ydata,lb,ub);
        p_win = msgbox('Parameter values for 3-rd degree c(T) estimated.', ...
            'User information','help');
        
        Max_err = max(abs(ydata-fun_deg_3(handles.p_est,tdata,handles)));
        Max_err_win = msgbox(['Maximal absoulte error |y-u|: ', ...
            num2str(Max_err)],'User information','help');
        
        set(handles.est_vals,'Visible','on')
        
        set(handles.p1_val,'Visible','on')
        set(handles.p2_val,'Visible','on')
        set(handles.p3_val,'Visible','on')
        set(handles.p4_val,'Visible','on')
        
        set(handles.p1_ed,'Visible','on')
        set(handles.p1_ed,'String',handles.p_est(1))
        set(handles.p2_ed,'Visible','on')
        set(handles.p2_ed,'String',handles.p_est(2))
        set(handles.p3_ed,'Visible','on')
        set(handles.p3_ed,'String',handles.p_est(3))
        set(handles.p4_ed,'Visible','on')
        set(handles.p4_ed,'String',handles.p_est(4))
        
        set(handles.export_temp,'Enable','on')
        set(handles.export_heat,'Enable','on')
        set(handles.plot_heat,'Enable','on')
        %r = fun_deg_3(p0,tdata,handles);
        %win = msgbox('The handles structure passed as an argument.', ...
        %                'User information','help');
        %r_win = msgbox('Numerical solution with 3-rd degree c(T) found.', ...
        %                'User information','help');
        %%%%%%%%%%%%%%%%%%%%%%%
        % Plot estimated c(T) %
        %%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.c_ax)
        handles.cplot = ...
            plot(handles.T-273.15,cp_eval(handles.p_est,handles.T));
        set(handles.cplot,'LineWidth',3,'LineStyle',':','Color','c')
        set(handles.c_ax,'FontSize',8)
        grid on
        xlabel('\bf{T}')
        ylabel('\bf{c(T)}')
        handles.hLeg_cp = legend('\it{3-rd degree c(T) approximation}');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Connecting with the context menu %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.cplot,'UIContextMenu',handles.line_properties)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot dimensional solution with % 
        % estimated parameters and the   %
        % experimental data  T(t)        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.ax_main)
        %cla(handles.ax_main)
        set(handles.ax_main,'FontSize',8)
        hold on
        grid on
        xlabel('\bf{t}')
        ylabel('\bf{T(t), [C deg]}')
        %%%%%%%%%%%%%%%%%%%%%
        % Experimental data %
        %%%%%%%%%%%%%%%%%%%%%
        handles.data_plot = ...
            plot(handles.t(handles.v), ...
            handles.T(handles.v)-273.15);
        set(handles.data_plot,'Color','b', ...
                              'LineStyle','-', ...
                              'LineWidth',3)
        handles.hLeg_data_plot = ...
            legend('\it{Experimental data}');
        %%%%%%%%%%%%%%%%%%%%%%
        % Numerical solution %
        %%%%%%%%%%%%%%%%%%%%%%
        handles.num_sol = ...
            fun_deg_3(handles.p_est,tdata,handles);
        %axes(handles.ax_main)
        %hold on
        handles.num_plot = ...
            plot(handles.t(handles.v),handles.num_sol-273.15);
        set(handles.num_plot,'Color','g', ...
                             'LineStyle',':', ...
                             'LineWidth',3)
        handles.hLeg_num_plot = ...
            legend('\it{Experimental data}','\it{Numerical solution}');
        %%%%%%%%%%%%%%%%%%%%%%
        %    Save handles    %
        %%%%%%%%%%%%%%%%%%%%%%
        guidata(hObject,handles);
    case 4
        axes(handles.ax_main)
        cla(handles.ax_main)
        
        axes(handles.c_ax)
        cla(handles.c_ax)
        %=============================================
        % input dialogue box, 4-th degree c(T)
        %=============================================
        %set(hObject,'Enable','off');
        %%%%%%%%%%%%%%%%
        %    Bounds    %
        %%%%%%%%%%%%%%%%
        bnd_prompt = {'Lower bounds for [p_1,p_2,p_3,p_4,p_5]', ...
          'Upper bounds for [p_1,p_2,p_3,p_4,p_5]'};
        bnd_wtitle = ...
            'Coefficient bounds for c(T) = p_{5}*T^4 + p_{4}*T^3 + p_{3}*T^2 + p_{2}*T + p_{1}';
        bnd_dims = [1,130];
        bnd_definput = {'',''};
        opts.Interpreter = 'tex';
        bnd = ...
            inputdlg(bnd_prompt,bnd_wtitle,bnd_dims,bnd_definput,opts);
        lb = str2num(bnd{1});
        handles.lb = lb;
        ub = str2num(bnd{2});
        handles.ub = ub;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check for correct input values %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while(numel(lb)~=5 || numel(ub)~=5)
            err = errordlg(['Out of bounds. Please try with correct ', ...
                'number of input values.'],'Input error');
            pause(3);
            %close(err);
            bnd_prompt = {'Lower bounds for [p_1,p_2,_3,p_4,p_5]', ...
                  'Upper bounds for [p_1,p_2,p_3,p_4,p_5]'};
            bnd_wtitle = ...
                'Coefficient bounds for c(T) = p_{5}*T^4 + p_{4}*T^3 + p_{3}*T^2 + p_{2}*T + p_{1}';
            bnd_dims = [1,130];
            bnd_definput = {'',''};
            opts.Interpreter = 'tex';
            bnd = inputdlg(bnd_prompt,bnd_wtitle,bnd_dims,bnd_definput,opts);
            lb = str2num(bnd{1});
            ub = str2num(bnd{2});
        end
        %%%%%%%%%%%%%%%%%%%%%
        % Extracting bounds %
        %%%%%%%%%%%%%%%%%%%%%
        lb = str2num(bnd{1});
        handles.lb = lb;
        ub = str2num(bnd{2});
        handles.ub = ub;
        %%%%%%%%%%%%%%%%
        % Coefficients %
        %%%%%%%%%%%%%%%%
        prompt = {'p_{1}','p_{2}','p_{3}','p_{4}','p_{5}'};
        wtitle = ...
            'Initial coefficients of c(T) = p_{5}*T^4 + p_{4}*T^3 + p_{3}*T^2 + p_{2}*T + p_{1}';
        dims = [1,125];
        definput = {'','','','',''};
        opts.Interpreter = 'tex';
        x = inputdlg(prompt,wtitle,dims,definput,opts);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Initial values for p0  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        p0(1) = str2num(x{1});
        p0(2) = str2num(x{2});
        p0(3) = str2num(x{3});
        p0(4) = str2num(x{4});
        p0(5) = str2num(x{5});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Estimating p vector    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %handles.fun = @(t,u) 2*u+3*t;
        %w_info = msgbox('Function f(t,u) successfully created.', ...
        %                'User information','help');
        %p = p0; %[1,1,1,1,1];
        %tdata = handles.t_dimless; %linspace(0,1,numel(handles.t));
        anon_f_4 = @(p,tdata) fun_deg_4(p,tdata,handles);
        tdata = handles.t_dimless; %linspace(0,1,numel(handles.t));
        ydata = handles.T(handles.start_i:handles.end_i)'; %fun_deg_4(p0,tdata,handles);
        handles.p_est = lsqcurvefit(anon_f_4,p0,tdata,ydata,lb,ub);
        p_win = msgbox('Parameter values for 4-th degree c(T) estimated.', ...
            'User information','help');
        
        Max_err = max(abs(ydata-fun_deg_4(handles.p_est,tdata,handles)));
        Max_err_win = msgbox(['Maximal absoulte error |y-u|: ', ...
            num2str(Max_err)],'User information','help');
        
        set(handles.est_vals,'Visible','on')
        
        set(handles.p1_val,'Visible','on')
        set(handles.p2_val,'Visible','on')
        set(handles.p3_val,'Visible','on')
        set(handles.p4_val,'Visible','on')
        set(handles.p5_val,'Visible','on')
        
        set(handles.p1_ed,'Visible','on')
        set(handles.p1_ed,'String',handles.p_est(1))
        set(handles.p2_ed,'Visible','on')
        set(handles.p2_ed,'String',handles.p_est(2))
        set(handles.p3_ed,'Visible','on')
        set(handles.p3_ed,'String',handles.p_est(3))
        set(handles.p4_ed,'Visible','on')
        set(handles.p4_ed,'String',handles.p_est(4))
        set(handles.p5_ed,'Visible','on')
        set(handles.p5_ed,'String',handles.p_est(5))
        
        set(handles.export_temp,'Enable','on')
        set(handles.export_heat,'Enable','on')
        set(handles.plot_heat,'Enable','on')
        %r = fun_deg_4(p0,tdata,handles);
        %win = msgbox('The handles structure passed as an argument.', ...
        %                'User information','help');
        %r_win = msgbox('Numerical solution with 4-th degree c(T) found.', ...
        %                'User information','help');
        %%%%%%%%%%%%%%%%%%%%%%%
        % Plot estimated c(T) %
        %%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.c_ax)
        handles.cplot = ...
            plot(handles.T-273.15,cp_eval(handles.p_est,handles.T));
        set(handles.cplot,'LineWidth',3,'LineStyle',':','Color','m')
        set(handles.c_ax,'FontSize',8)
        grid on
        xlabel('\bf{T}')
        ylabel('\bf{c(T)}')
        handles.hLeg_cp = legend('\it{4-th degree c(T) approximation}');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Connecting with the context menu %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.cplot,'UIContextMenu',handles.line_properties)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot dimensional solution with % 
        % estimated parameters and the   %
        % experimental data  T(t)        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.ax_main)
        cla(handles.ax_main)
        set(handles.ax_main,'FontSize',8)
        hold on
        grid on
        xlabel('\bf{t}')
        ylabel('\bf{T(t), [C deg]}')
        %%%%%%%%%%%%%%%%%%%%%
        % Experimental data %
        %%%%%%%%%%%%%%%%%%%%%
        handles.data_plot = ...
            plot(handles.t(handles.v), ...
            handles.T(handles.v)-273.15);
        set(handles.data_plot,'Color','b', ...
                              'LineStyle','-', ...
                              'LineWidth',3)
        handles.hLeg_data_plot = ...
            legend('\it{Experimental data}');
        %%%%%%%%%%%%%%%%%%%%%%
        % Numerical solution %
        %%%%%%%%%%%%%%%%%%%%%%
        handles.num_sol = ...
            fun_deg_4(handles.p_est,tdata,handles);
        %axes(handles.ax_main)
        %hold on
        handles.num_plot = ...
            plot(handles.t(handles.v),handles.num_sol-273.15);
        set(handles.num_plot,'Color','g', ...
                             'LineStyle',':', ...
                             'LineWidth',3)
        handles.hLeg_num_plot = ...
            legend('\it{Experimental data}','\it{Numerical solution}');
        %%%%%%%%%%%%%%%%%%%%%%
        %    Save handles    %
        %%%%%%%%%%%%%%%%%%%%%%
        guidata(hObject,handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIGHT-HAND SIDE FUNCTIONS FOR ESTIMATING THE VALUES  %
%          OF THE PARAMETERS IN THE FOUR CASES         %
%                FOR THE DEGREE OF c(T)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-st degree                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = fun_deg_1(p,tdata,handles)    
    % function f(t,u) for ode45
    fun = @(t,u) ...
	     (handles.tend * handles.eps * ...
 		  handles.S * handles.sigma * ...
         (M_fun(handles.M_args,handles.tend*t).^4 - u.^4) + ...
          handles.tend * handles.k * ...
		  handles.S * ...
		 (M_fun(handles.M_args,handles.tend*t) - u)) ./ ...
         ((p(1) + p(2)*u) * handles.m + ...
          c0_fun(handles.c0_args,u) * handles.m0 + ...
          handles.cs * handles.ms);
    
    [T,Y] = ode45(fun,tdata,handles.yinit); 
    % tdata - dimensionless!!!
    res = Y';%Y(:,1);
    %res = res'; % END OF THE FUNCTION
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-nd degree                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = fun_deg_2(p,tdata,handles)    
    % function f(t,u) for ode45
    fun = @(t,u) ...
         (handles.tend * handles.eps * ...
          handles.S * handles.sigma * ...
         (M_fun(handles.M_args,handles.tend*t).^4 - u.^4) + ...
          handles.tend * handles.k * ...
          handles.S * ...
         (M_fun(handles.M_args,handles.tend*t) - u)) ./ ...
         ((p(1) + p(2)*u + p(3)*u.^2) * handles.m + ...
          c0_fun(handles.c0_args,u) * handles.m0 + ...
          handles.cs * handles.ms);
    
    [T,Y] = ode45(fun,tdata,handles.yinit); 
    % tdata - dimensionless!!!
    res = Y';%Y(:,1);
    %res = res'; % END OF THE FUNCTION
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-rd degree                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = fun_deg_3(p,tdata,handles)    
    % function f(t,u) for ode45
    fun = @(t,u) ...
         (handles.tend * handles.eps * ...
          handles.S * handles.sigma * ...
         (M_fun(handles.M_args,handles.tend*t).^4 - u.^4) + ...
          handles.tend * handles.k * ...
          handles.S * ...
         (M_fun(handles.M_args,handles.tend*t) - u)) ./ ...
         ((p(1) + p(2)*u + p(3)*u.^2 + p(4)*u.^3) * handles.m + ...
          c0_fun(handles.c0_args,u) * handles.m0 + ...
          handles.cs * handles.ms);
    
    [T,Y] = ode45(fun,tdata,handles.yinit); 
    % tdata - dimensionless!!!
    res = Y';%Y(:,1);
    %res = res'; % END OF THE FUNCTION
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4-th degree                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = fun_deg_4(p,tdata,handles)    
    % function f(t,u) for ode45
    fun = @(t,u) ...
         (handles.tend * handles.eps * ...
          handles.S * handles.sigma * ...
         (M_fun(handles.M_args,handles.tend*t).^4 - u.^4) + ...
          handles.tend * handles.k * ...
          handles.S * ...
         (M_fun(handles.M_args,handles.tend*t) - u)) ./ ...
         ((p(1) + p(2)*u + p(3)*u.^2 + p(4)*u.^3 + p(5)*u.^4) * ...
          handles.m + ...
          c0_fun(handles.c0_args,u) * handles.m0 + ...
          handles.cs * handles.ms);
    
    [T,Y] = ode45(fun,tdata,handles.yinit); 
    % tdata - dimensionless!!!
    res = Y';%Y(:,1);
    %res = res'; % END OF THE FUNCTION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --- Executes during object creation, after setting all properties.
function cp_fun_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cp_fun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p1_ed_Callback(hObject, eventdata, handles)
% hObject    handle to p1_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p1_ed as text
%        str2double(get(hObject,'String')) returns contents of p1_ed as a double


% --- Executes during object creation, after setting all properties.
function p1_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p1_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p2_ed_Callback(hObject, eventdata, handles)
% hObject    handle to p2_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p2_ed as text
%        str2double(get(hObject,'String')) returns contents of p2_ed as a double


% --- Executes during object creation, after setting all properties.
function p2_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p2_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p3_ed_Callback(hObject, eventdata, handles)
% hObject    handle to p3_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p3_ed as text
%        str2double(get(hObject,'String')) returns contents of p3_ed as a double


% --- Executes during object creation, after setting all properties.
function p3_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p3_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p4_ed_Callback(hObject, eventdata, handles)
% hObject    handle to p4_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p4_ed as text
%        str2double(get(hObject,'String')) returns contents of p4_ed as a double


% --- Executes during object creation, after setting all properties.
function p4_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p4_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start_index_Callback(hObject, eventdata, handles)
% hObject    handle to start_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_index as text
%        str2double(get(hObject,'String')) returns contents of start_index as a double
tab_dat = get(handles.table_data,'Data');
sz_tab_dat = size(tab_dat);
ndata = sz_tab_dat(1);

start_i = get(hObject,'String');
%display(start_i)
%class(start_i)
if(isempty(start_i) || ...
       ~isnumeric(str2num(start_i)) || ...
       strcmp(num2str(str2num(start_i)),'') || ...
       str2num(start_i)<=0 || ...
       str2num(start_i)>ndata || ...
       floor(str2num(start_i))~=ceil(str2num(start_i)))
            f = errordlg(['Start time index is not valid.', ...
                ' Please try with a correct numerical value!'], ...
                'Input error');
end

if(~isempty(start_i) && ...
       isnumeric(str2num(start_i)) && ...
       ~strcmp(num2str(str2num(start_i)),'') && ...
       str2num(start_i)>0 && ... 
       str2num(start_i)<=ndata && ...
       floor(str2num(start_i))==ceil(str2num(start_i)))
            set(hObject,'Enable','off');
end

handles.start_i = str2num(start_i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Saving handles structure     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function start_index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function end_index_Callback(hObject, eventdata, handles)
% hObject    handle to end_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_index as text
%        str2double(get(hObject,'String')) returns contents of end_index as a double
tab_dat = get(handles.table_data,'Data');
sz_tab_dat = size(tab_dat);
ndata = sz_tab_dat(1);
%display(ndata)

end_i = get(hObject,'String');
%display(end_i)
%class(end_i)
if(isempty(end_i) || ...
       ~isnumeric(str2num(end_i)) || ...
       strcmp(num2str(str2num(end_i)),'') || ...
       str2num(end_i)<=0 || ...
       str2num(end_i)>ndata || ...
       floor(str2num(end_i))~=ceil(str2num(end_i)))
            f = errordlg(['End time index is not valid.', ...
                ' Please try with a correct numerical value!'], ...
                'Input error');
end

if(~isempty(end_i) && ...
       isnumeric(str2num(end_i)) && ...
       ~strcmp(num2str(str2num(end_i)),'') && ...
       str2num(end_i)>0 && ...
       str2num(end_i)<=ndata && ...
       floor(str2num(end_i))==ceil(str2num(end_i)))
            set(hObject,'Enable','off');
end

handles.end_i = str2num(end_i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Saving handles structure     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function end_index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p5_ed_Callback(hObject, eventdata, handles)
% hObject    handle to p5_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p5_ed as text
%        str2double(get(hObject,'String')) returns contents of p5_ed as a double


% --- Executes during object creation, after setting all properties.
function p5_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p5_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function line_color_Callback(hObject, eventdata, handles)
% hObject    handle to line_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function line_width_Callback(hObject, eventdata, handles)
% hObject    handle to line_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function line_style_Callback(hObject, eventdata, handles)
% hObject    handle to line_style (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function line_properties_Callback(hObject, eventdata, handles)
% hObject    handle to line_properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function solid_Callback(hObject, eventdata, handles)
% hObject    handle to solid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cplot,'LineStyle','-')


% --------------------------------------------------------------------
function dashed_Callback(hObject, eventdata, handles)
% hObject    handle to dashed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cplot,'LineStyle','--')

% --------------------------------------------------------------------
function dotted_Callback(hObject, eventdata, handles)
% hObject    handle to dotted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cplot,'LineStyle',':')

% --------------------------------------------------------------------
function one_width_Callback(hObject, eventdata, handles)
% hObject    handle to one_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cplot,'LineWidth',1.0)

% --------------------------------------------------------------------
function two_width_Callback(hObject, eventdata, handles)
% hObject    handle to two_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cplot,'LineWidth',2.0)

% --------------------------------------------------------------------
function three_width_Callback(hObject, eventdata, handles)
% hObject    handle to three_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cplot,'LineWidth',3.0)

% --------------------------------------------------------------------
function blue_color_Callback(hObject, eventdata, handles)
% hObject    handle to blue_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cplot,'Color','b')

% --------------------------------------------------------------------
function red_color_Callback(hObject, eventdata, handles)
% hObject    handle to red_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cplot,'Color','r')

% --------------------------------------------------------------------
function green_color_Callback(hObject, eventdata, handles)
% hObject    handle to green_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cplot,'Color','g')

% --------------------------------------------------------------------
function cyan_color_Callback(hObject, eventdata, handles)
% hObject    handle to cyan_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cplot,'Color','cyan')

% --------------------------------------------------------------------
function magenta_color_Callback(hObject, eventdata, handles)
% hObject    handle to magenta_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cplot,'Color','m')


% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function plot_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function plot_heat_Callback(hObject, eventdata, handles)
% hObject    handle to plot_heat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = questdlg('Is an analytic formula for heat capacity known?', ...
	'Heat capacity formula', ...
	'Yes','No','Yes');
if(strcmp(answer,'Yes'))
    % dlgbox for typing the function c(T)
    prompt = {'Enter an analytic heat capacity formula (T in K deg):'};
    dlgtitle = 'Input';
    dims = [1,50];
    res = inputdlg(prompt,dlgtitle,dims);
    fun = inline(res{1});
    % plot
    figure(3)
    %%%
    plot(handles.num_sol-273.15,fun(handles.num_sol-273.15),'b','LineWidth',3)
    hold on, grid on
    plot(handles.num_sol-273.15,cp_eval(handles.p_est,handles.num_sol),'r:','LineWidth',3)
    set(gca,'FontSize',16)
    xlabel('\bf{Temperature T, [C]}')
    ylabel('\bf{Heat capacity c(T), [J/(kg.C)]}')
    legend('\bf{Exact c(T)}','\bf{Estimated c(T)}')
    % info
    win = msgbox('Plotting heat capacity graph has been performed successfully.', ...
        'User information','help');
else
    win = msgbox(['Plotting analytic heat capacity graph is not possible since ', ...
        'an analytic formula has not been provided.'], ...
        'User information','error');
end


% --------------------------------------------------------------------
function export_temp_Callback(hObject, eventdata, handles)
% hObject    handle to export_temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A = [(handles.t(handles.v)); (handles.num_sol)];
%display(size(handles.t(handles.v)))
%display(size(handles.num_sol))
fid = fopen('solution.txt','w');
fprintf(fid,'%10f %10f\n',A);
fclose(fid);
mb = msgbox(['Numerical solution of the ODE corresponding to ', ...
             'the estimated heat coefficients values exported to ', ...
             'a text file named solution.txt.'],'User information','help');


% --------------------------------------------------------------------
function export_heat_Callback(hObject, eventdata, handles)
% hObject    handle to export_heat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cp_poly_val = get(handles.cp_fun,'Value');
fid = fopen('p_estimated.txt','w');
switch cp_poly_val
    case 1
        % 1-st degree polynomial
        fprintf(fid,'%10f %10f',handles.p_est);
        fclose(fid);
        mb = msgbox(['Estimated heat coefficients values exported to ', ...
            'a text file named p_estimated.txt.'],'User information','help');
    case 2
        % 2-nd degree polynomial
        fprintf(fid,'%10f %10f %10f',handles.p_est);
        fclose(fid);
        mb = msgbox(['Estimated heat coefficients values exported to ', ...
            'a text file named p_estimated.txt.'],'User information','help');
    case 3
        % 3-rd degree polynomial
        fprintf(fid,'%10f %10f %10f %10f',handles.p_est);
        fclose(fid);
        mb = msgbox(['Estimated heat coefficients values exported to ', ...
            'a text file named p_estimated.txt.'],'User information','help');
    case 4
        % 4-th degree polynomial
        fprintf(fid,'%10f %10f %10f %10f %10f',handles.p_est);
        fclose(fid);
        mb = msgbox(['Estimated heat coefficients values exported to ', ...
            'a text file named p_estimated.txt.'],'User information','help');
end
