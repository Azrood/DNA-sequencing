from collections import Counter

import plotly
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.http.response import HttpResponse
from django.shortcuts import render, redirect

from pages.dna_analysis import analysis, analysis_plot
from pages.forms import FilesUploadForm, SignupForm, SigninForm
from pages.models import FilesUpload, Utilisateur


# Create your views here.

@login_required(login_url="login")
def profile(request):
    return render(request, "profile.html")

@login_required(login_url="login")
def compare(request):
    return render(request, "base.html")

@login_required(login_url="login")
def home_view(request):
    return render(request,"index.html")

@login_required(login_url="login")
def analyse_view(request,*args,**kwargs):
    if request.method == 'POST':
        form = FilesUploadForm(request.POST,request.FILES)
        if form.is_valid():
            form.save()
            return render(request,"analyse.html",{'form': form})
    else:
        form = FilesUploadForm()
    print(Counter(analysis().seq))

    return render(request,"analyse.html",{'form': form,'dna_desc':analysis().description,'nuc_freq':dict(Counter(analysis().seq)),'graph_div':plotly.offline.plot(analysis_plot(), auto_open = False, output_type="div")})

@login_required(login_url="login")
def vis3d_view(request,*args,**kwargs):
    return render(request,"3dvis.html",{})

def register(request):
    if request.method == 'POST':
        form = SignupForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            password = form.cleaned_data.get('password1')
            user = authenticate(username=username, password=password)
            Utilisateur(user=user).save()
            login(request, user)
            return redirect('/home/')
    else:
        form = SignupForm()
    return render(request,'register.html/',{'form':form})

@login_required()
def disconnect(request):
    logout(request)
    return redirect('/login/')


def login_view(request):
    if request.user.is_authenticated:
        return redirect("/home/")
    if request.method == 'POST':
        form = SigninForm(data=request.POST)
        if form.is_valid():
            user = authenticate(
                    username=form.cleaned_data['username'],
                    password=form.cleaned_data['password'])
            login(request, user)
            return redirect('/home/')
    else:
        form = SigninForm()
    return render(request, "login.html", {'form': form})
