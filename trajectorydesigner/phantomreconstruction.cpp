#include "phantomreconstruction.h"

//#include "phantom.h"

PhantomReconstruction::PhantomReconstruction(QWidget *parent) : QWidget(parent)
{
	setAutoFillBackground(true);

	QPalette palette;
	palette.setColor(QPalette::Background, Qt::black);
	setPalette(palette);
}
