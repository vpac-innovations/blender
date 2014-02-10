/*
 * Copyright 2011, Blender Foundation.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Contributor: 
 *		Jeroen Bakker 
 *		Monique Dewanchand
 */

#include "COM_DifferenceMatteNode.h"
#include "BKE_node.h"
#include "COM_DifferenceMatteOperation.h"
#include "COM_SetAlphaOperation.h"

DifferenceMatteNode::DifferenceMatteNode(bNode *editorNode) : Node(editorNode)
{
	/* pass */
}

void DifferenceMatteNode::convertToOperations(NodeCompiler *compiler, const CompositorContext *context) const
{
	InputSocket *inputSocket = this->getInputSocket(0);
	InputSocket *inputSocket2 = this->getInputSocket(1);
	OutputSocket *outputSocketImage = this->getOutputSocket(0);
	OutputSocket *outputSocketMatte = this->getOutputSocket(1);
	bNode *editorNode = this->getbNode();

	DifferenceMatteOperation *operationSet = new DifferenceMatteOperation();
	operationSet->setSettings((NodeChroma *)editorNode->storage);
	compiler->addOperation(operationSet);
	
	compiler->mapInputSocket(inputSocket, operationSet->getInputSocket(0));
	compiler->mapInputSocket(inputSocket2, operationSet->getInputSocket(1));
	compiler->mapOutputSocket(outputSocketMatte, operationSet->getOutputSocket(0));

	SetAlphaOperation *operation = new SetAlphaOperation();
	compiler->addOperation(operation);
	
	compiler->mapInputSocket(inputSocket, operation->getInputSocket(0));
	compiler->addConnection(operationSet->getOutputSocket(), operation->getInputSocket(1));
	compiler->mapOutputSocket(outputSocketImage, operation->getOutputSocket());
	
	compiler->addOutputPreview(operation->getOutputSocket());
}
