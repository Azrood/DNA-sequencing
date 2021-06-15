from django import forms
from django.core.exceptions import ValidationError
from pages.models import FilesUpload, Utilisateur, User
from django.contrib.auth.forms import AuthenticationForm

CSS_CLASS_FIELD = "u-border-2 u-border-black u-border-no-left u-border-no-right u-border-no-top u-input u-input-rectangle u-white"

class FilesUploadForm(forms.ModelForm):
    class Meta:
        model = FilesUpload
        fields = ['file']

class SignupForm(forms.Form):
    username = forms.CharField(max_length=150, label='Nom d\'utilisateur', widget=forms.TextInput(attrs={'placeholder':'Nom d\'utilisateur', 'class':CSS_CLASS_FIELD}))
    password1 = forms.CharField(label='Mot de passe', widget=forms.PasswordInput(
            attrs={'placeholder':'Mot de passe',
            'pattern':r"(?=.*\d)(?=.*[a-z])(?=.*[A-Z]).{8,}", 
            'title':"Doit contenir au moins un nombre, une lettre majuscule et une lettre minuscule, et au moins 8 caractères",
            'class':CSS_CLASS_FIELD}))
    password2 = forms.CharField(label='Confirmer le mot de passe', widget=forms.PasswordInput(attrs={'placeholder':'Confirmez votre mot de passe', 'class':CSS_CLASS_FIELD}))
    email = forms.EmailField(label='Adresse Email',widget=forms.EmailInput(attrs={'placeholder':'Email', 'class':CSS_CLASS_FIELD}))

    def clean_username(self):
        username = self.cleaned_data['username'].lower()
        r = Utilisateur.user.get_queryset().filter(username=username)
        if r.count():
            raise ValidationError("Nom d'utilisateur déjà pris")
        return username
    
    def clean_email(self):
        email = self.cleaned_data['email'].lower()
        r = Utilisateur.user.get_queryset().filter(email=email)
        if r.count():
            raise ValidationError("Cette adresse existe déjà")
        return email

    def clean_password2(self):
        password1 = self.cleaned_data.get('password1')
        password2 = self.cleaned_data.get('password2')
        if password1 and password2 and password1 != password2:
            raise forms.ValidationError("Les mots de passes ne correspondent pas")
        return password2

    def save(self, commit=True):
        user = User.objects.create_user(
            self.cleaned_data['username'],
            self.cleaned_data['email'],
            self.cleaned_data['password1']
        )
        if commit:
            user.save()
        return user

class SigninForm(AuthenticationForm):
    username = forms.CharField(max_length=150, label='Nom d\'utilisateur', widget=forms.TextInput(attrs={'placeholder':'Entrez le nom d\'utilisateur', 'class':CSS_CLASS_FIELD}))
    password = forms.CharField(label='Mot de passe', widget=forms.PasswordInput(attrs={'placeholder':'Entrez votre mot de passe', 'class':CSS_CLASS_FIELD}))

    error_messages = {
        'invalid_login': ("Nom d`'utilisateur ou mot de passe incorrect"),
    }