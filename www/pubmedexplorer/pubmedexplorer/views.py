from django.shortcuts import render
from django.http import HttpResponse

def index(request):
    return render(
        request,
        "pubmedexplorer/index.html"
    )

def search(request):
    return HttpResponse("Hello, world. You're at the pubmedexplorer search page.")