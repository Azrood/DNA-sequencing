# Generated by Django 3.2.4 on 2021-06-15 00:42

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pages', '0005_alter_filesupload_file'),
    ]

    operations = [
        migrations.AlterField(
            model_name='filesupload',
            name='file',
            field=models.FileField(upload_to='media/'),
        ),
    ]
