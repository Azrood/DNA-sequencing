from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from django.forms import ModelForm


# Create your models here.
class FilesUpload(models.Model):
    file = models.FileField(upload_to='media/')
    type_fichier = models.TextField("Type de fichier", max_length=30, blank=True)

class Utilisateur(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE)
    fichiers = models.ForeignKey(FilesUpload, on_delete=models.SET_NULL, blank=True, null=True)
