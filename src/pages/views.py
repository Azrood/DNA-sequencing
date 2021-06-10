from django.http.response import HttpResponse
from django.shortcuts import render
from django.http import HttpResponse
from pages.dna_analysis import analysis
from pages.models import FilesUpload
from .forms import FilesUploadFrom
from .dna_analysis import analysis
# Create your views here.
def home_view(request,*args,**kwargs):
    return render(request,"index.html",{})

def analyse_view(request,*args,**kwargs):
    if request.method == 'POST':
        form = FilesUploadFrom(request.POST,request.FILES)
        if form.is_valid():
            form.save()
            return render(request,"analyse.html",{'form': form})
    else:
        form = FilesUploadFrom()

    return render(request,"analyse.html",{'form': form,'dna_desc':analysis()})

def vis3d_view(request,*args,**kwargs):
    return render(request,"3dvis.html",{})