from django import forms
from django.core.exceptions import ValidationError
from pages.models import FilesUpload, Utilisateur, User
from django.contrib.auth.forms import AuthenticationForm, PasswordChangeForm

CSS_CLASS_LOGIN_FIELD = "u-border-2 u-border-black u-border-no-left u-border-no-right u-border-no-top u-input u-input-rectangle u-white"
CSS_FORM_PASS_MAIL_CHANGE = "u-border-1 u-border-grey-30 u-input u-input-rectangle u-white"

class FilesUploadForm(forms.ModelForm):
    type_fichier = forms.CharField(max_length=30, label="type_fichier", widget=forms.TextInput(attrs={"placeholder": "Type du fichier"}))
    class Meta:
        model = FilesUpload
        fields = ['file', 'type_fichier']
    
    def save(self, user, commit=True):
        form = FilesUpload.objects.create(
            file=self.cleaned_data['file'],
            type_fichier=self.cleaned_data['type_fichier'],
            utilisateur=Utilisateur.objects.get(user=user)
        )
        if commit:
            form.save()
        return form

class SignupForm(forms.Form):
    username = forms.CharField(max_length=150, label='Nom d\'utilisateur', widget=forms.TextInput(attrs={'placeholder':'Nom d\'utilisateur', 'class':CSS_CLASS_LOGIN_FIELD}))
    password1 = forms.CharField(label='Mot de passe', widget=forms.PasswordInput(
            attrs={'placeholder':'Mot de passe',
            'pattern':r"(?=.*\d)(?=.*[a-z])(?=.*[A-Z]).{8,}", 
            'title':"Doit contenir au moins un nombre, une lettre majuscule et une lettre minuscule, et au moins 8 caractères",
            'class':CSS_CLASS_LOGIN_FIELD}))
    password2 = forms.CharField(label='Confirmer le mot de passe', widget=forms.PasswordInput(attrs={'placeholder':'Confirmez votre mot de passe', 'class':CSS_CLASS_LOGIN_FIELD}))
    email = forms.EmailField(label='Adresse Email',widget=forms.EmailInput(attrs={'placeholder':'Email', 'class':CSS_CLASS_LOGIN_FIELD}))

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
    username = forms.CharField(max_length=150, label='Nom d\'utilisateur', widget=forms.TextInput(attrs={'placeholder':'Entrez le nom d\'utilisateur', 'class':CSS_CLASS_LOGIN_FIELD}))
    password = forms.CharField(label='Mot de passe', widget=forms.PasswordInput(attrs={'placeholder':'Entrez votre mot de passe', 'class':CSS_CLASS_LOGIN_FIELD}))

    error_messages = {
        'invalid_login': ("Nom d'utilisateur ou mot de passe incorrect"),
    }

class ChangePass(PasswordChangeForm):
    new_password1 = forms.CharField(
        label="Nouveau mot de passe",
        widget=forms.PasswordInput(attrs={'autocomplete': 'new-password','placeholder':'Entrez votre nouveau mot de passe', 'class':CSS_FORM_PASS_MAIL_CHANGE}),
        strip=False,
        help_text="Au minimum 8 caractères.\nDoit contenir au moins un caractère numérique, une majuscule et un symbole",
    )
    new_password2 = forms.CharField(
        label="Confirmation du nouveau mot de passe",
        strip=False,
        widget=forms.PasswordInput(attrs={'autocomplete': 'new-password','placeholder':'Confirmer votre nouveau mot de passe', 'class':CSS_FORM_PASS_MAIL_CHANGE}),
    )
    old_password = forms.CharField(
        label="Mot de passe actuel",
        strip=False,
        widget=forms.PasswordInput(attrs={'autocomplete': 'current-password', 'autofocus': True,'placeholder':'Entrer votre mot de passe actuel', 'class':CSS_FORM_PASS_MAIL_CHANGE})
    )

    error_messages = {
        'password_mismatch':('Les mots de passe ne correspondent pas'),
        'password_incorrect':("Votre ancien mot de passe est incorrect"),
    }

class ChangeMail(forms.Form):
    error_messages={
        'email_inuse': "Cette adresse est déjà enregistrée.",
        'password_incorrect': 'Mot de passe incorrect'
}

    old_mail = forms.EmailField(label='Adresse Email actuelle',widget=forms.EmailInput(attrs={'placeholder':'Entrer votre adresse actuelle', 'class':CSS_FORM_PASS_MAIL_CHANGE, "disabled":True}))
    new_email = forms.EmailField(label='Nouvelle adresse Email',widget=forms.EmailInput(attrs={'placeholder':'Entrez votre nouvelle adresse email', 'class':CSS_FORM_PASS_MAIL_CHANGE}))
    current_password = forms.CharField(label="Mot de passe", widget=forms.PasswordInput(attrs={'placeholder':'entrez votre mot de passe pour confirmer', 'class':CSS_FORM_PASS_MAIL_CHANGE}),required=True)

    def __init__(self, user, *args, **kwargs):
        self.user = user
        super(ChangeMail, self).__init__(*args, **kwargs)
    
    def clean_current_password(self):
        """
        Validates that the password field is correct.
        """
        current_password = self.cleaned_data["current_password"]
        if not self.user.check_password(current_password):
            raise forms.ValidationError(self.error_messages['password_incorrect'], code='password_incorrect',)
        return current_password
    
    def clean_new_email(self):
        """
        Prevents an e-mail address that is already registered from being registered by a different user.
        """
        email = self.cleaned_data.get('new_email')
        if Utilisateur.user.get_queryset().filter(email=email).count() > 0:
            raise forms.ValidationError(self.error_messages['email_inuse'], code='email_inuse',)
        return email

    def clean_old_mail(self):
        return self.cleaned_data['old_mail']

    def save(self, commit=True):
        self.user.email=self.cleaned_data["new_email"]
        self.user.save()
        return self.user