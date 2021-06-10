from django import forms
from .models import FilesUpload

class FilesUploadFrom(forms.ModelForm):
    class Meta:
        model = FilesUpload
        fields = ['file']