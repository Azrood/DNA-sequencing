from django.db import models
from django.forms import ModelForm
from django.conf import settings
# Create your models here.
class FilesUpload(models.Model):
    file = models.FileField(upload_to='media/')


