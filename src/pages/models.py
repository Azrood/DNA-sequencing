from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from django.dispatch import receiver
from django.forms import ModelForm
from django.db.models.signals import post_delete

# Create your models here.

class Utilisateur(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE)

class FilesUpload(models.Model):
    file = models.FileField(upload_to="media/")
    type_fichier = models.TextField("Type de fichier", max_length=30, blank=True)
    utilisateur = models.ForeignKey(Utilisateur, on_delete=models.SET_NULL, blank=True, null=True)

@receiver(post_delete, sender=FilesUpload)
def mymodel_delete(sender, instance, **kwargs):
    # Pass false so FileField doesn't save the model.
    instance.file.delete(False)
