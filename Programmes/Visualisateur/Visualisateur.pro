#-------------------------------------------------
#
# Project created by QtCreator 2017-04-30T00:18:51
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Visualisateur
TEMPLATE = app

INCLUDEPATH += /usr/include/qwt
LIBS += -lqwt-qt5

CONFIG += c++14


SOURCES += main.cpp\
        mainwindow.cpp

HEADERS  += mainwindow.h \
    wavelets.h

FORMS    += mainwindow.ui
