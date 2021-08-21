function varargout = ode_gui_solve(varargin)
% ODE_GUI_SOLVE MATLAB code for ode_gui_solve.fig
%      ODE_GUI_SOLVE, by itself, creates a new ODE_GUI_SOLVE or raises the existing
%      singleton*.
%
%      H = ODE_GUI_SOLVE returns the handle to a new ODE_GUI_SOLVE or the handle to
%      the existing singleton*.
%
%      ODE_GUI_SOLVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ODE_GUI_SOLVE.M with the given input arguments.
%
%      ODE_GUI_SOLVE('Property','Value',...) creates a new ODE_GUI_SOLVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ode_gui_solve_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ode_gui_solve_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ode_gui_solve

% Last Modified by GUIDE v2.5 23-Jun-2021 23:22:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ode_gui_solve_OpeningFcn, ...
                   'gui_OutputFcn',  @ode_gui_solve_OutputFcn, ...
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


% --- Executes just before ode_gui_solve is made visible.
function ode_gui_solve_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ode_gui_solve (see VARARGIN)

% Choose default command line output for ode_gui_solve
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ode_gui_solve wait for user response (see UIRESUME)
% uiwait(handles.figure1);
    
    %=============================================
    % title
    %=============================================
    axes(handles.stat_f)
    pos = get(handles.stat_f,'position');
    text(pos(3)/2-pos(3)/3,pos(4)/2,'$$ u'' = f(x,u), \, \, \, u(x_0)=u_0, \, \, \, x \in [x_0; \, X] $$', ...
        'interpreter','latex',...
        'FontSize',14);
    set(handles.stat_f,'TickDir','out');
    set(handles.stat_f,'visible','off');
    axis([0,pos(3),0,pos(4)])
    %=============================================
    % a
    %=============================================
    axes(handles.aAx)
    pos = get(handles.aAx,'position');
    text(pos(3)/2,pos(4)/2,'$$ x_0: $$','interpreter','latex',...
        'FontSize',14);
    set(handles.aAx,'TickDir','out');
    set(handles.aAx,'visible','off');
    axis([0,pos(3),0,pos(4)])
    %=============================================
    % b
    %=============================================
    axes(handles.bAx)
    pos = get(handles.bAx,'position');
    text(pos(3)/2,pos(4)/2,'$$ X: $$','interpreter','latex',...
        'FontSize',14);
    set(handles.bAx,'TickDir','out');
    set(handles.bAx,'visible','off');
    axis([0,pos(3),0,pos(4)])
    %=============================================
    % u0
    %=============================================
    axes(handles.uAx)
    pos = get(handles.uAx,'position');
    text(pos(3)/2,pos(4)/2,'$$ u_{0}: $$','interpreter','latex',...
        'FontSize',14);
    set(handles.uAx,'TickDir','out');
    set(handles.uAx,'visible','off');
    axis([0,pos(3),0,pos(4)])
    %=============================================
    % f(x,u)
    %=============================================
    axes(handles.fAx)
    pos = get(handles.fAx,'position');
    text(pos(3)/2,pos(4)/2,'$$ f(x,u): $$','interpreter','latex',...
        'FontSize',14);
    set(handles.fAx,'TickDir','out');
    set(handles.fAx,'visible','off');
    axis([0,pos(3),0,pos(4)])
    %=============================================
    % u(x)
    %=============================================
    axes(handles.usolAx)
    pos = get(handles.usolAx,'position');
    text(pos(3)/2,pos(4)/2,'$$ u(x): $$','interpreter','latex',...
        'FontSize',14);
    set(handles.usolAx,'TickDir','out');
    set(handles.usolAx,'visible','off');
    axis([0,pos(3),0,pos(4)])


% --- Outputs from this function are returned to the command line.
function varargout = ode_gui_solve_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function aVal_Callback(hObject, eventdata, handles)
% hObject    handle to aVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of aVal as text
%        str2double(get(hObject,'String')) returns contents of aVal as a double
as = get(handles.aVal,'String'); a = str2num(as);
bs = get(handles.bVal,'String'); b = str2num(bs);
u0s = get(handles.initVal,'String'); u0 = str2num(u0s);
r_side = get(handles.fExpr,'String'); f = inline(r_side,'x','u');

if(~isempty(as) && ~isempty(a) && ...
       ~isempty(bs) && ~isempty(b) && ...
       ~isempty(u0s) && ~isempty(u0) && ...
       ~isempty(r_side))
   set(handles.solveBut,'Enable','on');
end


% --- Executes during object creation, after setting all properties.
function aVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bVal_Callback(hObject, eventdata, handles)
% hObject    handle to bVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bVal as text
%        str2double(get(hObject,'String')) returns contents of bVal as a double
as = get(handles.aVal,'String'); a = str2num(as);
bs = get(handles.bVal,'String'); b = str2num(bs);
u0s = get(handles.initVal,'String'); u0 = str2num(u0s);
r_side = get(handles.fExpr,'String'); f = inline(r_side,'x','u');

if(~isempty(as) && ~isempty(a) && ...
       ~isempty(bs) && ~isempty(b) && ...
       ~isempty(u0s) && ~isempty(u0) && ...
       ~isempty(r_side))
   set(handles.solveBut,'Enable','on');
end


% --- Executes during object creation, after setting all properties.
function bVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initVal_Callback(hObject, eventdata, handles)
% hObject    handle to initVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initVal as text
%        str2double(get(hObject,'String')) returns contents of initVal as a double
as = get(handles.aVal,'String'); a = str2num(as);
bs = get(handles.bVal,'String'); b = str2num(bs);
u0s = get(handles.initVal,'String'); u0 = str2num(u0s);
r_side = get(handles.fExpr,'String'); f = inline(r_side,'x','u');

if(~isempty(as) && ~isempty(a) && ...
       ~isempty(bs) && ~isempty(b) && ...
       ~isempty(u0s) && ~isempty(u0) && ...
       ~isempty(r_side))
   set(handles.solveBut,'Enable','on');
end


% --- Executes during object creation, after setting all properties.
function initVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fExpr_Callback(hObject, eventdata, handles)
% hObject    handle to fExpr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fExpr as text
%        str2double(get(hObject,'String')) returns contents of fExpr as a double
as = get(handles.aVal,'String'); a = str2num(as);
bs = get(handles.bVal,'String'); b = str2num(bs);
u0s = get(handles.initVal,'String'); u0 = str2num(u0s);
r_side = get(handles.fExpr,'String'); f = inline(r_side,'x','u');

if(~isempty(as) && ~isempty(a) && ...
       ~isempty(bs) && ~isempty(b) && ...
       ~isempty(u0s) && ~isempty(u0) && ...
       ~isempty(r_side))
   set(handles.solveBut,'Enable','on');
end


% --- Executes during object creation, after setting all properties.
function fExpr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fExpr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in solveBut.
function solveBut_Callback(hObject, eventdata, handles)
% hObject    handle to solveBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

as = get(handles.aVal,'String'); a = str2num(as);
bs = get(handles.bVal,'String'); b = str2num(bs);
if(a>b)
    t = a;
    a = b;
    b = t;
end
u0s = get(handles.initVal,'String'); u0 = str2num(u0s);
r_side = get(handles.fExpr,'String'); f = inline(r_side,'x','u');

set(hObject,'Enable','off');

% Grid
x = linspace(a,b,300);
N = numel(x);
% Grid step
h = x(2)-x(1);
% Numeric solution
y = zeros(1,N);
% Initial condition
y(1) = u0;
% Improved Euler method
for i = 1 : N-1
    ym = y(i)+h*f(x(i),y(i));
         y(i+1) = y(i) + h/2 * ...
         (f(x(i),y(i)) + f(x(i+1),ym));
end
% Save in handles structure
handles.x = x; 
handles.y = y;
handles.h = h;
% Plot the results
axes(handles.axMain) % make current axes
Line = plot(handles.x,handles.y);
     set(Line,'Color','b');
     set(Line,'LineStyle','-');
     set(Line,'LineWidth',3);
handles.Line = Line;
hold on, grid on
guidata(gcbo,handles);

set(handles.clearBut,'Enable','on');
set(handles.anal_sol_script,'Enable','on');
set(handles.save_menu,'Enable','on')


% --- Executes on button press in clearBut.
function clearBut_Callback(hObject, eventdata, handles)
% hObject    handle to clearBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','off');
set(handles.aVal,'String','');
set(handles.bVal,'String','');
set(handles.initVal,'String','');
set(handles.fExpr,'String','');
%%%
res = get(handles.anal_sol,'String');
if(~isempty(res))
    set(handles.anal_sol,'String','');
    set(handles.anal_sol,'Enable','off');
    set(handles.anal_sol_script,'Enable','off');
    cla(handles.ax1,'reset'); grid off;
    cla(handles.ax2,'reset'); grid off;
    set(handles.uip,'Visible','off');
end
%%%
set(handles.anal_sol_script,'Enable','off');
set(handles.anal_sol,'Enable','off');
set(handles.err,'Enable','off');
cla(handles.axMain,'reset'); grid off;
%set(handles.save_menu,'Enable','off');



% --- Executes on button press in anal_sol_script.
function anal_sol_script_Callback(hObject, eventdata, handles)
% hObject    handle to anal_sol_script (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.anal_sol,'Enable','on');
ue = get(handles.anal_sol,'String');
if(~isempty(ue))
    u_fun = inline(ue,'x');
    u = u_fun(handles.x);
    handles.u = u;
    set(handles.err,'Enable','on');
    set(handles.anal_sol_script,'Enable','off');
    cla(handles.axMain,'reset'); grid off;
    set(handles.axMain,'Visible','off');
    %%%
    uip = uipanel('Position',[0.045,0.218,0.492,0.714]);
    set(uip,'BackgroundColor',[1,1,1]);
    set(uip,'FontSize',12);
    handles.uip = uip;
    %uip = uipanel('Position',[-0.037,0.122,0.634,0.876]);
    %%%
    ax1 = subplot(1,2,1,'Parent',uip);
    handles.ax1 = ax1;
    Line_u = plot(handles.x,handles.u,'b','LineWidth',3);
    handles.Line_u = Line_u;
    hold on, grid on;
    Line_y = plot(handles.x,handles.y,'c--','LineWidth',3);
    handles.Line_y = Line_y;
    legend('\it{\bf{Точно решение}}','\it{\bf{Приближено решение}}');
    %set(ax1,'FontSize',12);
    guidata(gcbo,handles);
    %%%
    ax2 = subplot(1,2,2,'Parent',uip);
    handles.ax2 = ax2;
    sol_err = abs(handles.u-handles.y);
    handles.sol_err = sol_err;
    Line_err = plot(handles.x,handles.sol_err,'r','LineWidth',3);
    hold on, grid on;
    legend('\it{\bf{Грешка}}')
    %set(ax2,'FontSize',12);
    handles.Line_err = Line_err;
    guidata(gcbo,handles);
end


function anal_sol_Callback(hObject, eventdata, handles)
% hObject    handle to anal_sol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of anal_sol as text
%        str2double(get(hObject,'String')) returns contents of anal_sol as a double


% --- Executes during object creation, after setting all properties.
function anal_sol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anal_sol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in err.
function err_Callback(hObject, eventdata, handles)
% hObject    handle to err (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns err contents as cell array
%        contents{get(hObject,'Value')} returns selected item from err
Num = get(hObject,'Value');
switch Num
    case 1 % Max. absolute error
        m_err = max(handles.sol_err);
        handles.m_err = m_err;
        helpdlg(['Максимална абсолютна грешка за подобрения метод на Ойлер: ', ...
            num2str(handles.m_err)], ...
            'Информация');
    case 2 % Relative error
        rel_err = max(handles.sol_err) / max(abs(handles.u));
        handles.rel_err = rel_err;
        helpdlg(['Относителна грешка за подобрения метод на Ойлер: ', ...
            num2str(handles.rel_err)], ...
            'Информация');
end


% --- Executes during object creation, after setting all properties.
function err_CreateFcn(hObject, eventdata, handles)
% hObject    handle to err (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function save_menu_Callback(hObject, eventdata, handles)
% hObject    handle to save_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xx = handles.x; yy = handles.y; hh = handles.h;
guidata(hObject,handles);
save('data_x.mat','xx');
save('data_y.mat','yy');
save('data_h.mat','hh');
set(hObject,'Enable','off');
