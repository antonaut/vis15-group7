/****************************************************************************
** Meta object code from reading C++ file 'ClassSelectionDialog.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../ClassSelectionDialog.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ClassSelectionDialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_ClassSelectionDialog[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      22,   21,   21,   21, 0x0a,
      54,   21,   21,   21, 0x0a,
      84,   72,   21,   21, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_ClassSelectionDialog[] = {
    "ClassSelectionDialog\0\0"
    "sltDirectTypeInChanged(QString)\0"
    "sltEnterPressed()\0item,column\0"
    "sltDoubleClicked(QTreeWidgetItem*,int)\0"
};

void ClassSelectionDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        ClassSelectionDialog *_t = static_cast<ClassSelectionDialog *>(_o);
        switch (_id) {
        case 0: _t->sltDirectTypeInChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: _t->sltEnterPressed(); break;
        case 2: _t->sltDoubleClicked((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData ClassSelectionDialog::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject ClassSelectionDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_ClassSelectionDialog,
      qt_meta_data_ClassSelectionDialog, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ClassSelectionDialog::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ClassSelectionDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ClassSelectionDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ClassSelectionDialog))
        return static_cast<void*>(const_cast< ClassSelectionDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int ClassSelectionDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
