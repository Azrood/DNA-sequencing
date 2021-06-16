from collections import Counter

import plotly
from django.contrib import messages
from django.contrib.auth import (authenticate, login, logout,
                                 update_session_auth_hash)
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.http.response import HttpResponse
from django.shortcuts import redirect, render

from pages.dna_analysis import analysis, analysis_plot, dotplotx
from pages.forms import (ChangeMail, ChangePass, FilesUploadForm, SigninForm,
                         SignupForm)
from pages.models import FilesUpload, Utilisateur
from Bio.SeqUtils import GC
import pandas as pd

# Create your views here.

@login_required(login_url="login")
def compare(request):
    return render(request, "compare.html",{'dotplot':dotplotx(analysis().seq[0:10],analysis().seq[0:10])})

@login_required(login_url="login")
def home_view(request):
    return render(request,"index.html")

@login_required(login_url="login")
def analyse_view(request,*args,**kwargs):
    utilisateur = Utilisateur.objects.get(user=request.user)
    fichiers = utilisateur.filesupload_set.all()
    if request.method == 'POST':
        form = FilesUploadForm(request.POST,request.FILES)
        file_fasta = FilesUpload.objects.get(pk=request.POST['files'])
        record = analysis(f"media/{file_fasta.file.name}")
        context = {'fichiers': fichiers, 'form':form,'dna':record,'nuc_freq':dict(Counter(record.seq)),'nuc_graph':plotly.offline.plot(analysis_plot(record.seq), auto_open = False, output_type="div"),'compo_adn':{'GC':GC(record.seq),'AT':100-GC(record.seq)},'transcribe':record.seq.transcribe(),'translate':record.seq.translate(),'complement':record.seq.complement(),'aa_freq':dict(Counter(record.seq.translate())),'aa_graph':plotly.offline.plot(analysis_plot(record.seq.translate()), auto_open = False, output_type="div")}
        return render(request, "analyse.html", context)
    else:
        form = FilesUploadForm()
    # amino_acids = [a for a in analysis().seq.translate().split('*')]
    # df = pd.DataFrame({'amino_acids':amino_acids})
    context = {"fichiers": fichiers}
    return render(request, "analyse.html", context)

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


@login_required(login_url="login")
def profile(request):
    utilisateur = Utilisateur.objects.get(user=request.user)
    fichiers = utilisateur.filesupload_set.all()
    
    file_form = FilesUploadForm()
    pass_form = ChangePass(user=request.user)
    mail_form = ChangeMail(user=request.user, initial={'old_mail':request.user.email})
    if request.method == 'POST':
        if "add" in request.POST:
            file_form = FilesUploadForm(data=request.POST, files=request.FILES)
            if file_form.is_valid():
                file_form.save(user=request.user)
        elif "del" in request.POST:
            ids =[x for x in request.POST.keys() if x not in ("csrfmiddlewaretoken", "del")]
            files = FilesUpload.objects.filter(id__in=ids)
            for f in files:
                f.delete()

        elif "pass_change" in request.POST:
            pass_form = ChangePass(user=request.user, data=request.POST)
            if pass_form.is_valid():
                pass_form.save()
                update_session_auth_hash(request, pass_form.user)
                messages.success(request, 'Votre mot de passe a été changé avec succès')
                redirect("/home/")
        elif "mail_change" in request.POST:
            mail_form = ChangeMail(user=request.user, data=request.POST)
            if mail_form.is_valid():
                mail_form.save()
                messages.success(request, 'Votre adresse mail a été modifiée avec succès')
                redirect("/home/")
    return render(request, 'profile.html', {'passform':pass_form,
                                            'mailform':mail_form,
                                            'fileform': file_form,
                                            'files': fichiers})
