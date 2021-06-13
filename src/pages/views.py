from django.http.response import HttpResponse
from django.shortcuts import render
from django.http import HttpResponse
from pages.dna_analysis import analysis
from pages.models import FilesUpload
from .forms import FilesUploadFrom
from .dna_analysis import analysis,analysis_plot
from collections import Counter
import plotly
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
    print(Counter(analysis().seq))

    return render(request,"analyse.html",{'form': form,'dna_desc':analysis().description,'nuc_freq':dict(Counter(analysis().seq)),'graph_div':plotly.offline.plot(analysis_plot(), auto_open = False, output_type="div")})

def vis3d_view(request,*args,**kwargs):
    return render(request,"3dvis.html",{})